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

#include "svenssonHaggkvistCanopySource.H"
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
    defineTypeNameAndDebug(svenssonHaggkvistCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        svenssonHaggkvistCanopySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::svenssonHaggkvistCanopySource::svenssonHaggkvistCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
  canopySource(name, modelType, dict, mesh),
  CpEps1_
  (
        dimensionedScalar::lookupOrAddToDict
        (
            "CpEps1",
            coeffs_,
            1.8
        )
  )
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::fv::svenssonHaggkvistCanopySource::addSup
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


void Foam::fv::svenssonHaggkvistCanopySource::addSup
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


void Foam::fv::svenssonHaggkvistCanopySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k")) {

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    eqn += canopy*pow(mag(U),3);
  }
  else if (eqn.psi().name() == word("epsilon")) {

    const volScalarField& epsilon = eqn.psi();
    const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    
    fvMatrix<scalar> Sepsilon
    (
       fvm::Sp(canopy*CpEps1_/k*pow(mag(U),3), epsilon)
    );

    eqn += Sepsilon;
  }
}


void Foam::fv::svenssonHaggkvistCanopySource::addSup
 (
     const volScalarField& rho,
     fvMatrix<scalar>& eqn,
     const label fieldi
 )
 {

  const volScalarField& canopy = canopy_;

  if (eqn.psi().name() == word("k")) {   
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    eqn += rho*canopy*pow(mag(U),3);
  }
  else if (eqn.psi().name() == word("epsilon")) {

    const volScalarField& epsilon = eqn.psi();
    const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");


    fvMatrix<scalar> Sepsilon
    (
     fvm::Sp(rho*canopy*CpEps1_/k*pow(mag(U),3), epsilon)
    );

    
    eqn += Sepsilon;    
  }
}



// ************************************************************************* //
