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
    Foam::fv::svenssonHaggkvistCanopySource

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
        type            svenssonHaggkvistCanopySource;
        active          on;
    }
    \endverbatim

SourceFiles
    svenssonHaggkvistCanopySource.C

\*---------------------------------------------------------------------------*/

#ifndef svenssonHaggkvistCanopySource_H
#define svenssonHaggkvistCanopySource_H

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
                 Class svenssonHaggkvistCanopySource Declaration
\*---------------------------------------------------------------------------*/

class svenssonHaggkvistCanopySource
:
    public canopySource
{
private:

        // Canopy coefficients
        dimensionedScalar CpEps1_;
  
        //- Disallow default bitwise copy construct
        svenssonHaggkvistCanopySource(const svenssonHaggkvistCanopySource&);

        //- Disallow default bitwise assignment
        void operator=(const svenssonHaggkvistCanopySource&);


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
    TypeName("svenssonHaggkvistCanopySource");


    // Constructors

        //- Construct from components
        svenssonHaggkvistCanopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~svenssonHaggkvistCanopySource()
    {}


    // Member Functions
        // access

       dimensionedScalar CpEps1() const
        {
            return CpEps1_;
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
