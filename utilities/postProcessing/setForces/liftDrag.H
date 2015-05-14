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

Description

\*---------------------------------------------------------------------------*/

#include "boundBox.H"

scalarField boxMask(const boundBox& bb, const pointField& points)
{
    scalarField mask(points.size(),1.0);
    if (bb.mag() > 0.0)
    {
        forAll (mask, i)
        {
            mask[i] = scalar(bb.contains(points[i]));
        }
    }
    return mask;
}


vector pressureLoad
(
    const volScalarField& p,
    const word& patchName,
    const boundBox bbox=boundBox()
)
{
    const fvMesh& mesh = p.mesh();

    label patchLabel = mesh.boundaryMesh().findPatchID(patchName);

    vector pressureForce;

    if (patchLabel != -1)
    {
        const pointField &centres=mesh.C().boundaryField()[patchLabel];
        scalarField mask = boxMask(bbox,centres);
        pressureForce = sum
        (
            mask*p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
    }

    return pressureForce;
}

vector pressureMoment
(
    const volScalarField& p,
    const word& patchName,
    const vector cPoint,
    const boundBox bbox=boundBox()
)
{
    const fvMesh& mesh = p.mesh();

    label patchLabel = mesh.boundaryMesh().findPatchID(patchName);

    vector pressureForce;
    vector moments;

    if (patchLabel != -1)
    {
        const pointField &centres=mesh.C().boundaryField()[patchLabel];
        scalarField mask = boxMask(bbox,centres);
        vectorField pLoads=mask*p.boundaryField()[patchLabel]*mesh.Sf().boundaryField()[patchLabel];
        vectorField R = centres-cPoint;
        moments = sum(pLoads^R);
    }
    else
    {
        moments = sum(vectorField(0));
    }

    return moments;
}


vector viscousLoad
(
    const volVectorField& U,
    const dimensionedScalar& mu,
    const word& patchName
)
{
    const fvMesh& mesh = U.mesh();

    label patchLabel = mesh.boundaryMesh().findPatchID(patchName);

    vector viscousForce;

    if (patchLabel != -1)
    {
        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );
    }
    else
    {
        viscousForce = sum(vectorField(0));
    }

    return viscousForce;
}


vector turbulentLoad
(
    const autoPtr<turbulenceModel>& turbulence,
    const volVectorField& U,
    const dimensionedScalar& mu,
    const word& patchName
)
{
    const fvMesh& mesh = U.mesh();

    label patchLabel = mesh.boundaryMesh().findPatchID(patchName);

    vector turbForce;

    if (patchLabel != -1)
    {
        turbForce = sum
        (
          - mesh.Sf().boundaryField()[patchLabel]
          & turbulence->R()().boundaryField()[patchLabel]
        );
    }
    else
    {
        turbForce = sum(vectorField(0));
    }

    return turbForce;
}

/*

vector momentCoefficient
(
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref,
    const scalar& Lref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "vector momentCoefficient\n"
            "(\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref,\n"
            "    const scalar& Lref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForceMoment;
    vector viscousForceMoment;

    if (patchLabel != -1)
    {
        pressureForceMoment = sum
        (
            mesh.Cf().boundaryField()[patchLabel]^
            (p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel])
        );

        viscousForceMoment = sum
        (
            mesh.Cf().boundaryField()[patchLabel]^
            (-mu.value()*U.boundaryField()[patchLabel].snGrad()
            *mesh.magSf().boundaryField()[patchLabel])
        );
    }
    else
    {
        pressureForceMoment = sum(vectorField(0));
        viscousForceMoment = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;
    scalar Mref = Fref*Lref;

    vector pressureCoeff = pressureForceMoment/Mref;
    vector viscousCoeff = viscousForceMoment/Mref;

    return (pressureCoeff + viscousCoeff);
}

*/

// ************************************************************************* //
