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
def create_fields( runTime, mesh ):
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info() << "Reading field p\n" << nl
    
    from Foam.OpenFOAM import IOobject, word, fileName
    from Foam.finiteVolume import volScalarField
    p = volScalarField( IOobject( word( "p" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    ext_Info() << "Reading field U\n" << nl
    from Foam.finiteVolume import volVectorField
    U = volVectorField( IOobject( word( "U" ),
                                  fileName( runTime.timeName() ),
                                  mesh,
                                  IOobject.MUST_READ,
                                  IOobject.AUTO_WRITE ),
                        mesh )

    from Foam.finiteVolume.cfdTools.incompressible import createPhi
    phi = createPhi( runTime, mesh, U )

    pRefCell = 0
    pRefValue = 0.0
    
    from Foam.finiteVolume import setRefCell
    pRefCell, pRefValue = setRefCell( p, mesh.solutionDict().subDict( word( "SIMPLE" ) ), pRefCell, pRefValue )

    from Foam.transportModels import singlePhaseTransportModel
    laminarTransport = singlePhaseTransportModel( U, phi )

    from Foam import incompressible
    turbulence = incompressible.RASModel.New( U, phi, laminarTransport )

    from Foam.finiteVolume import porousZones
    pZones = porousZones( mesh )
    
    from Foam.OpenFOAM import Switch
    pressureImplicitPorosity = Switch( False )

    nUCorr = 0
    if pZones.size():
       # nUCorrectors for pressureImplicitPorosity
       if (mesh.solutionDict().subDict( word( "SIMPLE" ) ).found( word( "nUCorrectors" ) ) ) :
          from Foam.OpenFOAM import readInt
          nUCorr = readInt( mesh.solutionDict().subDict( word( "SIMPLE" ) ).lookup( word( "nUCorrectors" ) ) )
          pass
       if nUCorr > 0 :
          pressureImplicitPorosity = True
          pass
       pass

    return p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence, pZones, pressureImplicitPorosity, nUCorr


#---------------------------------------------------------------------------
def initConvergenceCheck( simple ):
    eqnResidual = 1
    maxResidual = 0
    convergenceCriterion = 0.0
    
    from Foam.OpenFOAM import word
    tmp, convergenceCriterion = simple.readIfPresent( word( "convergence" ), convergenceCriterion )

    return eqnResidual, maxResidual, convergenceCriterion


#---------------------------------------------------------------------------
def convergenceCheck( runTime, maxResidual, convergenceCriterion ):
    if maxResidual < convergenceCriterion :
        from Foam.OpenFOAM import ext_Info, nl
        ext_Info << "reached convergence criterion: " << convergenceCriterion <<nl
        runTime.writeAndEnd()
        ext_Info << "latestTime =" << runTime.timeName()
        pass
    
    pass
    

#-----------------------------------------------------------------------------------------
def fun_UEqn( mesh, phi, U, p, turbulence, pZones, nUCorr, pressureImplicitPorosity, eqnResidual, maxResidual ):
    
    from Foam import fvm, fvc    
    # Construct the Momentum equation

    # The initial C++ expression does not work properly, because of
    #  1. turbulence.divDevRhoReff( U ) - changes values for the U boundaries
    #  2. the order of expression arguments computation differs with C++
    #UEqn = fvm.div( phi, U ) + turbulence.divDevReff( U ) 

    UEqn = turbulence.divDevReff( U ) + fvm.div( phi, U ) 

    UEqn.relax()

    # Include the porous media resistance and solve the momentum equation
    # either implicit in the tensorial resistance or transport using by
    # including the spherical part of the resistance in the momentum diagonal

    trAU = None
    trTU = None
    if pressureImplicitPorosity :
        from Foam.OpenFOAM import I, tensor
        tTU = tensor( I ) * UEqn.A()

        pZones.addResistance( UEqn, tTU )
    
        trTU = tTU.inv()
               
        from Foam.OpenFOAM import word
        trTU.rename( word( "rAU" ) )
        
        for UCorr in range ( nUCorr ):
            U.ext_assign( trTU & ( UEqn.H() - fvc.grad( p ) ) )
            pass
        
        U.correctBoundaryConditions()
        pass
    else:
        pZones.addResistance( UEqn )
        
        from Foam.finiteVolume import solve
        eqnResidual =  solve( UEqn == - fvc.grad( p ) ).initialResidual()
        
        maxResidual = max( eqnResidual, maxResidual )
        
        trAU = 1.0 / UEqn.A()
        
        from Foam.OpenFOAM import word
        trAU.rename( word( "rAU" ) )
        pass
    
    return UEqn, trTU, trAU, eqnResidual, maxResidual


#---------------------------------------------------------------------------
def fun_pEqn( mesh, p, U, trTU, trAU, UEqn, phi, runTime, pressureImplicitPorosity, nNonOrthCorr, \
              eqnResidual, maxResidual, cumulativeContErr, pRefCell, pRefValue, ):
   
    if pressureImplicitPorosity :
       U.ext_assign( trTU & UEqn.H() )
       pass
    else:
       U.ext_assign( trAU * UEqn.H() )
       pass
    
    UEqn.clear() 
    
    from Foam import fvc, fvm
    phi.ext_assign( fvc.interpolate( U ) & mesh.Sf() )

    from Foam.finiteVolume import adjustPhi
    adjustPhi( phi, U, p )
    
    for nonOrth in range( nNonOrthCorr + 1 ) :
        tpEqn = None
        if pressureImplicitPorosity :
            tpEqn = ( fvm.laplacian( trTU, p ) == fvc.div( phi ) )
            pass
        else:
            tpEqn = ( fvm.laplacian( trAU, p ) == fvc.div( phi ) )
            pass
        
        tpEqn.setReference( pRefCell, pRefValue )
        # retain the residual from the first iteration
        if nonOrth == 0 :
            eqnResidual = tpEqn.solve().initialResidual()
            maxResidual = max( eqnResidual, maxResidual )
            pass
        else:
            tpEqn.solve()
            pass
        
        if nonOrth == nNonOrthCorr :
            phi.ext_assign( phi - tpEqn.flux() )
            pass
        
        pass
    
    from Foam.finiteVolume.cfdTools.incompressible import continuityErrs
    cumulativeContErr = continuityErrs( mesh, phi, runTime, cumulativeContErr )

    # Explicitly relax pressure for momentum corrector
    p.relax()
           
    if pressureImplicitPorosity :
        U.ext_assign( U - ( trTU & fvc.grad( p ) ) )
    else:
        U.ext_assign( U - ( trAU * fvc.grad( p ) ) )
        pass
       
    U.correctBoundaryConditions()

    return eqnResidual, maxResidual

#---------------------------------------------------------------------------
def main_standalone( argc, argv ):

    from Foam.OpenFOAM.include import setRootCase
    args = setRootCase( argc, argv )

    from Foam.OpenFOAM.include import createTime
    runTime = createTime( args )

    from Foam.OpenFOAM.include import createMesh
    mesh = createMesh( runTime )

    p, U, phi, pRefCell, pRefValue, laminarTransport, turbulence, pZones, pressureImplicitPorosity, nUCorr = create_fields( runTime, mesh )
    
    from Foam.finiteVolume.cfdTools.general.include import initContinuityErrs
    cumulativeContErr = initContinuityErrs()
    
    from Foam.OpenFOAM import ext_Info, nl
    ext_Info()<< "\nStarting time loop\n" << nl
    
    while runTime.loop() :
        
        ext_Info() << "Time = " << runTime.timeName() << nl << nl
        
        from Foam.finiteVolume.cfdTools.general.include import readSIMPLEControls
        simple, nNonOrthCorr, momentumPredictor, transonic = readSIMPLEControls( mesh )
        
        eqnResidual, maxResidual, convergenceCriterion = initConvergenceCheck( simple )

        p.storePrevIter()


        UEqn, trTU, trAU, eqnResidual, maxResidual = fun_UEqn( mesh, phi, U, p, turbulence, pZones, nUCorr, pressureImplicitPorosity, eqnResidual, maxResidual )

        eqnResidual, maxResidual = fun_pEqn( mesh, p, U, trTU, trAU, UEqn, phi, runTime, pressureImplicitPorosity, nNonOrthCorr, \
                                             eqnResidual, maxResidual, cumulativeContErr, pRefCell, pRefValue, )

        turbulence.correct()

        runTime.write()
        
        ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << nl
        
        convergenceCheck( runTime, maxResidual, convergenceCriterion)
        pass
    
    ext_Info() << "End\n"

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "010701" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam1.7.1 or higher \n "


#--------------------------------------------------------------------------------------
