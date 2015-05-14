/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    speciesFoam

Description
    Steady-state solver for incompressible, turbulent 
    dispersion of a free number of species
Note
    Code based on the standard solver simpleFoam and pollutionFoam by 
    Albert Passalacqua (alberto.passalacqua@tin.it) 
    Molecular dispersion neglected.
Author
    David Segersson, SMHI(Swedish Meteorological And Hydrological Institute)
    david.segersson@smhi.se
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

    simpleControl simple(mesh);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //Setting schemes for species fractions
    // div discretisation scheme
    word divScheme ("div(phi,Yi)");

    // Laplacian discretisation scheme
    word laplacianScheme ("laplacian(diff,Yi)");
    
    Info<< "\nStarting time loop\n" << endl;

    //phi = fvc::interpolate(U) & mesh.Sf();
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	Info<<"Solving transport equation"<<endl;

	// Solving pollutant specie transport equation
	for(label i = 0; i < Y.size(); i++)
	  {      
	    volScalarField& Yi = Y[i];
	    volScalarField& eSTi = eST[i];

	    //	    volScalarField kappaEffConc
	    //  (
	    //   "kappaEffConc",
	    //    turbulence->nu()/PrConc + turbulence->nut()/PrConct
	    //    );
	    
		
	    fvScalarMatrix concEqn
	      (
	       fvm::div(phi, Yi, divScheme)
	       - fvm::laplacian((turbulence->nut())/SC_T, Yi, laplacianScheme)
	       ==
	       eSTi/densityAir
	       );

	    concEqn.relax();
	    concEqn.solve();
	    
	    //	    Yi.max(0.0);
	  }
        
        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    
    //Calculating max value in domain
    for(label i = 0; i < Y.size(); i++)
      {      
	volScalarField& Yi = Y[i];
	scalar maxMassFraction=0;
	
	forAll(Yi.internalField() ,cell_i)
	  {
	    if(Yi.internalField()[cell_i]>maxMassFraction)
	      maxMassFraction=Yi.internalField()[cell_i];
	  }   
		
	Info<< "Max mass fraction for "<<speciesNames[i]<<": "<<maxMassFraction<< endl;
	Info<< "Max concentration for "<<speciesNames[i]<<": "<<maxMassFraction*densityAir*1e9<<"microg/m3"<<endl;
      }

    Info<< "End\n" << endl;	
    
    return(0);
}


// ************************************************************************* //
