/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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

    Copyright is held by the author.

Class
    Foam::fv::dalpeMassonCanopySource

Group
    grpFvOptionsSources

Author
    David Segersson, Swedish Meteorological and Hydrological Institute, david.segersson@smhi.se

Description
    Momentum and turbulence source/sink-terms for tree canopy according to
    B. Dalpé and C. Masson. Numerical simulation of wind flow near a forest
    edge. J. of Wind Eng. and Ind. Aerodyn, 2006

Dalpé & Masson, 2006
    
Usage
    Example usage:
    \verbatim
    SBM
    {
        type            dalpeMassonCanopySource;
        active          on;
    }
    \endverbatim

SourceFiles
    dalpeMassonCanopySource.C

\*---------------------------------------------------------------------------*/

#ifndef dalpeMassonCanopySource_H
#define dalpeMassonCanopySource_H

#include "fvOption.H"
#include "dimensionedTypes.H"
#include "wallFvPatch.H"
#include "canopySource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class dalpeMassonCanopySource Declaration
\*---------------------------------------------------------------------------*/

class dalpeMassonCanopySource
:
    public canopySource
{
private:

        // Canopy coefficients
        dimensionedScalar betaP_;
        dimensionedScalar betaD_;
        dimensionedScalar C4_;
        dimensionedScalar C5_;
  
        //- Disallow default bitwise copy construct
        dalpeMassonCanopySource(const dalpeMassonCanopySource&);

        //- Disallow default bitwise assignment
        void operator=(const dalpeMassonCanopySource&);


        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        // Private Member Functions

  
public:


    //- Runtime type information
    TypeName("dalpeMassonCanopySource");


    // Constructors

        //- Construct from components
        dalpeMassonCanopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~dalpeMassonCanopySource()
    {}


    // Member Functions
        // access

       dimensionedScalar betaP() const
        {
            return betaP_;
        }

       dimensionedScalar betaD() const
        {
            return betaD_;
        }

       dimensionedScalar C4() const
        {
            return C4_;
        }

       dimensionedScalar C5() const
        {
            return C5_;
        }


        //- Source term to momentum equation
        void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to scalar equations (k & epsilon)
        void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        //- Source term to compressible momentum equation
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to compressible scalar equation (k & epsilon)
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
