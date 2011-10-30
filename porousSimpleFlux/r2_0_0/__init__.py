#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------
from Foam import ref, man

#---------------------------------------------------------------------------
def create_fields( runTime, mesh ):
    ref.ext_Info() << "Reading field p\n" << ref.nl
    
    p = man.volScalarField( man.IOobject( ref.word( "p" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    phi = man.createPhi( runTime, mesh, U )

    pRefCell = 0
    pRefValue = 0.0
    
    pRefCell, pRefValue = ref.setRefCell( p, mesh.solutionDict().subDict( ref.word( "SIMPLE" ) ), pRefCell, pRefValue )

    laminarTransport = man.singlePhaseTransportModel( U, phi )

    turbulence = man.incompressible.RASModel.New( U, phi, laminarTransport )


    return p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence 

#---------------------------------------------------------------------------
def createPorousZones( mesh, simple ):

     pZones = man.porousZones( mesh )
     pressureImplicitPorosity = ref.Switch( False )

     # nUCorrectors used for pressureImplicitPorosity
     nUCorr = 0
     
     if pZones.size():
         # nUCorrectors for pressureImplicitPorosity
         nUCorr = simple.dict().lookupOrDefault( ref.word( "nUCorrectors" ), 0 )

         if nUCorr > 0:
             pressureImplicitPorosity = ref.Switch( True );
             ref.ext_Info() << "Using pressure implicit porosity" << ref.nl
             pass
         else:
             ref.ext_Info() << "Using pressure explicit porosity" << ref.nl
             pass
     
     return pZones, pressureImplicitPorosity, nUCorr
    
    
#---------------------------------------------------------------------------
def fun_UEqn( mesh, phi, U, p, turbulence, pZones, nUCorr, pressureImplicitPorosity ):
    
    # Construct the Momentum equation

    # The initial C++ expression does not work properly, because of
    #  1. turbulence.divDevRhoReff( U ) - changes values for the U boundaries
    #  2. the order of expression arguments computation differs with C++
    #UEqn = fvm.div( phi, U ) + turbulence.divDevReff( U ) 

    UEqn = man.fvVectorMatrix( turbulence.divDevReff( U ), man.Deps( turbulence, U ) ) + man.fvm.div( phi, U ) 

    UEqn.relax()

    # Include the porous media resistance and solve the momentum equation
    # either implicit in the tensorial resistance or transport using by
    # including the spherical part of the resistance in the momentum diagonal

    trAU = None
    trTU = None
    if pressureImplicitPorosity :
        tTU = man.volTensorField( ref.tensor( ref.I ) * UEqn.A(), man.Deps( UEqn ) )

        pZones.addResistance( UEqn, tTU )
    
        trTU = man.volTensorField( tTU.inv(), man.Deps( tTU ) )
               
        trTU.rename( ref.word( "rAU" ) )
        
        for UCorr in range ( nUCorr ):
            U << ( trTU() & ( UEqn.H() - ref.fvc.grad( p ) ) ) # mixed calculations
            pass
        
        U.correctBoundaryConditions()
        pass
    else:
        pZones.addResistance( UEqn )
        
        ref.solve( UEqn == -man.fvc.grad( p ) )
        
        trAU = man.volScalarField( 1.0 / UEqn.A(), man.Deps( UEqn ) )
        
        trAU.rename( ref.word( "rAU" ) )
        pass
    
    return UEqn, trTU, trAU
    

#---------------------------------------------------------------------------
def fun_pEqn( mesh, simple, p, U, trTU, trAU, UEqn, phi, runTime, pressureImplicitPorosity, \
              cumulativeContErr, pRefCell, pRefValue, ):
   
    if pressureImplicitPorosity :
       U << ( trTU() & UEqn.H() ) # mixed calculations
       pass
    else:
       U <<  trAU * UEqn.H()
       pass
    
    #UEqn.clear() 
    
    phi << ( ref.fvc.interpolate( U ) & mesh.Sf() )

    ref.adjustPhi( phi, U, p )
    
    for nonOrth in range( simple.nNonOrthCorr() + 1 ) :
        tpEqn = None
        if pressureImplicitPorosity :
            tpEqn = ( ref.fvm.laplacian( trTU, p ) == ref.fvc.div( phi ) )
            pass
        else:
            tpEqn = ( ref.fvm.laplacian( trAU, p ) == ref.fvc.div( phi ) )
            pass
        
        tpEqn.setReference( pRefCell, pRefValue )
        # retain the residual from the first iteration
        if nonOrth == 0 :
            tpEqn.solve()
            pass
        else:
            tpEqn.solve()
            pass
        
        if nonOrth == simple.nNonOrthCorr() :
            phi-= tpEqn.flux()
            pass
        
        pass
    
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )

    # Explicitly relax pressure for momentum corrector
    p.relax()
           
    if pressureImplicitPorosity :
        U -= ( trTU() & ref.fvc.grad( p ) ) # mixed calcaulations
    else:
        U -= trAU * ref.fvc.grad( p )
        pass
       
    U.correctBoundaryConditions()

    return cumulativeContErr

#---------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    simple = man.simpleControl( mesh )

    p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence = create_fields( runTime, mesh )
    
    pZones, pressureImplicitPorosity, nUCorr = createPorousZones( mesh, simple )
    
    cumulativeContErr = ref.initContinuityErrs()
    
    ref.ext_Info()<< "\nStarting time loop\n" << ref.nl
    
    while simple.loop() :
        
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        
        p.storePrevIter()

        UEqn, trTU, trAU = fun_UEqn( mesh, phi, U, p, turbulence, pZones, nUCorr, pressureImplicitPorosity )

        cumulativeContErr = fun_pEqn( mesh, simple, p, U, trTU, trAU, UEqn, phi, runTime, pressureImplicitPorosity, \
                                             cumulativeContErr, pRefCell, pRefValue, )

        turbulence.correct()

        runTime.write()
        
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass
    
    ref.ext_Info() << "End\n"

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.0.0 or higher \n "


#--------------------------------------------------------------------------------------
