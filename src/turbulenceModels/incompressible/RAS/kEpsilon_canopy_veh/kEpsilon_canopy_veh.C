/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kEpsilon_canopy_veh.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilon_canopy_veh, 0);
addToRunTimeSelectionTable(RASModel, kEpsilon_canopy_veh, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilon_canopy_veh::kEpsilon_canopy_veh
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cfcar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cfcar",
            coeffDict_,
            0.64 //Mochida 2008
        )
    ),
    L_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "L",
            coeffDict_,
            10, //Mochida 2008
	    dimensionSet(0,1,0,0,0,0,0)			 
        )
    ),
    Cecar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cecar",
            coeffDict_,
            2.47  //Mochida 2008
        )
    ),
    betaP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaP",
            coeffDict_,
            1.0 //Dalpe and Masson 2007 
        )
    ),
    betaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaD",
            coeffDict_,
            5.03 //Dalpe and Masson 2007 
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            0.78//Dalpe and Masson 2007   
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            coeffDict_,
            0.78 //Dalpe and Masson 2007 
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
     ),
    canopy_
    (
        IOobject
        (
            "canopy",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
         mesh_
    ),
    Ucar_
    (
        IOobject
        (
            "Ucar",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    VAD_
    (
        IOobject
        (
            "VAD",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )


{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kEpsilon_canopy_veh::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kEpsilon_canopy_veh::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kEpsilon_canopy_veh::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> kEpsilon_canopy_veh::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool kEpsilon_canopy_veh::read()
{
    if (RASModel::read())
    {
        Cfcar_.readIfPresent(coeffDict());
        Cecar_.readIfPresent(coeffDict());
        L_.readIfPresent(coeffDict());
        betaP_.readIfPresent(coeffDict_);
        betaD_.readIfPresent(coeffDict_);
        C4_.readIfPresent(coeffDict_);
        C5_.readIfPresent(coeffDict_);
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilon_canopy_veh::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));

    //Source terms for canopy layer
    volScalarField Fk1("RASModel::Fk1",canopy_*betaP_*pow(mag(U_),3));
    volScalarField Fk2("RASModel::Fk2",canopy_*betaD_*k_*mag(U_));

    //Source terms for vehicle canopy
    volVectorField Fi = 0.5*Cfcar_*VAD_*mag(U_-Ucar_)*(U_-Ucar_);
    volScalarField Fk = (U_-Ucar_)&Fi;
    volScalarField Fepsilon = epsilon_*sqrt(k_)*Cecar_/L_;


    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Fepsilon
        + (C1_*G + C4_*Fk1)*epsilon_/k_
	- fvm::Sp((C5_*Fk2 + C2_*epsilon_)/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G+Fk1+Fk
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
