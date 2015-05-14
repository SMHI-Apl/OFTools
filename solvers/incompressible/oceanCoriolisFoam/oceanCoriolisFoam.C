/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    oceanCoriolisFoam

Description
    Transient solver for incompressible, turbulent flow of seawater with weakly 
    varying density (Boussinesq approximation). Originally based on turbFoam, 
    with added equations for temperature and salinity transport. Source terms 
    in the momentum equation for buoyancy and coriolis forces are included.
    
Author
    Walter Gyllenram

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RASModel/RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

#   include "myReadTransportProperties.H"
#   include "myReadEnvironmentalProperties.H"
#   include "myReadDensityCoeffs.H"
#   include "myReadLowerBoundValues.H"

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
 	     == (rho-a0)*g/a0 - (twoOmega ^ U)
           );

            if (momentumPredictor)
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rUA = 1.0/UEqn.A();

                U = rUA*UEqn.H();
                phi = (fvc::interpolate(U) & mesh.Sf()) 
                    + fvc::ddtPhiCorr(rUA, U, phi);

                adjustPhi(phi, U, p);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
		     fvm::laplacian(rUA, p) == fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

	// Solve energy equation
	tmp<volScalarField> myGammaEffT
	  (
	   nu/PrT + turbulence->nut()/PrTt
	   );
	
	tmp<fvScalarMatrix> TEqn
	  (
	   fvm::ddt(T)
	   + fvm::div(phi, T)
	   - fvm::laplacian(myGammaEffT, T)
	   );
	
	myGammaEffT.clear();
	
	TEqn().relax();
	
	solve(TEqn);
	
	bound(T,T0);
	
	// Solve salt transport equation
	tmp<volScalarField> myGammaEffS
	  (
	   nu/PrS + turbulence->nut()/PrSt
	   );
	
	tmp<fvScalarMatrix> SEqn
	  (
	   fvm::ddt(S)
	   + fvm::div(phi, S)
	   - fvm::laplacian(myGammaEffS, S)
	   );

	myGammaEffS.clear();

	SEqn().relax();
	
	solve(SEqn);
	
	bound(S,S0);

	volScalarField TC=T-T0C;	
	
	// Equation of state (UNESCO)
	/*
	rho = (a0 + a1*(TC) + a2*pow((TC),2) + a3*pow((TC),3) + a4*pow((TC),4) + a5*pow((TC),5))
	  + (b0 + b1*(TC) + b2*pow((TC),2) + b3*pow((TC),3) + b4*pow((TC),4))*S 
	  + (c0 + c1*(TC) + c2*pow((TC),2))*pow(S,3/2)
	  + d0*pow(S,2);
	*/

	rho = a0 + (a1 + (a2 + (a3 + (a4 + a5*TC)*TC)*TC)*TC)*TC
	  + (b0 + (b1 + (b2 + (b3 + b4*TC)*TC)*TC)*TC)*S
	  + (c0 + (c1 + c2*TC)*TC)*pow(S,3/2) + d0*S*S;
	
	TC.clear();

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
