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

\*---------------------------------------------------------------------------*/

#include "groundDist.H"
#include "patchWave.H"
#include "fvMesh.H"
#include "wallPolyPatch.H"
#include "fvPatchField.H"
#include "Field.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::groundDist::groundDist(const fvMesh& mesh, labelHashSet patchIDs, const bool correctWalls)
:
    volScalarField
    (
        IOobject
        (
            "y",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    cellDistFuncs(mesh),
    correctWalls_(correctWalls),
    nUnset_(0)
{
    groundDist::correct(patchIDs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::groundDist::~groundDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct for mesh geom/topo changes. Might be more intelligent in the
// future (if only small topology change)
void Foam::groundDist::correct(labelHashSet patchIDs)
{
    // Get patchids of walls
    //labelHashSet wallPatchIDs(getPatchIDs(wallPolyPatch::typeName));

    // Calculate distance starting from wallPatch faces.
    patchWave wave(cellDistFuncs::mesh(), patchIDs, correctWalls_);

    // Transfer cell values from wave into *this 
    transfer(wave.distance());

    // Transfer values on patches into boundaryField of *this
    forAll(boundaryField(), patchI)
    {
        if (boundaryField()[patchI].type() != emptyFvPatchScalarField::typeName)
        {
            scalarField& waveFld = wave.patchDistance()[patchI];

            boundaryField()[patchI].transfer(waveFld);
        }
    }

    // Transfer number of unset values
    nUnset_ = wave.nUnset();
}


// ************************************************************************* //
