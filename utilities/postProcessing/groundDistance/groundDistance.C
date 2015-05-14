/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  Application
  groundDistance

  Description
  Calculates a scalar field for distance to given patches

  \*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "fvCFD.H"
#include "groundDist.H"

#ifndef namespaceFoam
#define namespaceFoam    
using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readLanduseDict.H"

  labelHashSet patchIDs(sourcePatches.size());
      
    // Finding indexes to given patches 
    forAll(sourcePatches,sp)
      {
	label patchI;
	patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
	
	Info<< "Found source face "<< sourcePatches[sp]
	    <<" at patch with index " << patchI << endl;  
	
	if (patchI == -1)
	  {
	    FatalErrorIn(args.executable())
	      << "Cannot find patch "<< sourcePatches[sp] << exit(FatalError);
	  }

	patchIDs.insert(patchI);
      }

    volScalarField d
    (
        IOobject
        (
            "d",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("d", dimLength, SMALL),
        calculatedFvPatchScalarField::typeName
    );

    d = (const volScalarField&) groundDist(mesh,patchIDs).y();
    
    d.write();

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
