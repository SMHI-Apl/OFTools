/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "dalpeMassonCanopySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"
#include "fvOption.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(dalpeMassonCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        dalpeMassonCanopySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::dalpeMassonCanopySource::dalpeMassonCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
  canopySource(name, modelType, dict, mesh),
  betaP_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaP",
            coeffs_,
            1.0
        )
  ),
  betaD_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaD",
            coeffs_,
            5.03
        )
  ),
  C4_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "C4",
            coeffs_,
            0.78
        )
  ),
  C5_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "C5",
            coeffs_,
            0.78
        )
   )
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::fv::dalpeMassonCanopySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
  
    const volVectorField& U = eqn.psi();
    const volScalarField& canopy = canopy_;

    fvMatrix<vector> S_canopy
    (
       fvm::Sp(canopy * mag(U), U)
     );
    
    eqn -=  S_canopy;
  
}


void Foam::fv::dalpeMassonCanopySource::addSup
 (
     const volScalarField& rho,
     fvMatrix<vector>& eqn,
     const label fieldi
 )
 {

  if (eqn.psi().name() == word("U")) {
    
    const volVectorField& U = eqn.psi();
    const volScalarField& canopy = canopy_;

    fvMatrix<vector> S_canopy
    (
       fvm::Sp(rho*canopy * mag(U), U)
    );

    eqn -=  S_canopy;
  }

  
 }


void Foam::fv::dalpeMassonCanopySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k")) {

    // volScalarField Fk1("Fk1",canopy_*betaP_*pow(mag(U_),3));
    
    const volScalarField& k = eqn.psi();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    fvMatrix<scalar> Sk
    (
     betaP_*canopy*pow(mag(U),3) - fvm::Sp(betaD_*canopy*mag(U), k)
    );

    eqn +=  Sk;
  }
  else if (eqn.psi().name() == word("epsilon")) {

    const volScalarField& epsilon = eqn.psi();
    const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    
    fvMatrix<scalar> Sepsilon
    (
     fvm::Sp(canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), epsilon)
    );

    eqn += Sepsilon;    
  }
}


void Foam::fv::dalpeMassonCanopySource::addSup
 (
     const volScalarField& rho,
     fvMatrix<scalar>& eqn,
     const label fieldi
 )
 {

  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k")) {

    // volScalarField Fk1("Fk1",canopy_*betaP_*pow(mag(U_),3));
    
    const volScalarField& k = eqn.psi();
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    fvMatrix<scalar> Sk
    (
     betaP_*rho*canopy*pow(mag(U),3) - fvm::Sp(betaD_*rho*canopy*mag(U), k)
    );

    eqn +=  Sk;
  }
  else if (eqn.psi().name() == word("epsilon")) {

    const volScalarField& epsilon = eqn.psi();
    const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");


    fvMatrix<scalar> Sepsilon
    (
     fvm::SuSp(rho*canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), epsilon)
    );

    
    eqn += Sepsilon;    
  }
}



// ************************************************************************* //
