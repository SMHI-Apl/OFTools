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
    David Segersson, Swedish Meteorological and Hydrological Institute
    david.segersson@smhi.se

Description
    Momentum and turbulence source/sink-terms for tree canopy according to Svensson & Häggkvist (1990), with parameters adjusted according t    o Mochida et al. (2008).
    
    Momentum source term is expressed as:
        \f[
            S_u = -\rho C_d \alpha u^2
        \f]

    Source term for k is expressed as:

        \f[
            S_k = \rho C_d \alpha u^3
        \f]

    Source term for epsilon is expressed as:

        \f[
            S_{\epsilon} = \rho C_d \alpha C_{p \epsilon 1} \frac{\epsilon}{k} u^3 
        \f]

    where
    
    \vartable
        \rho             | density
        u                | magnitude of velocity
        k                | turbulent kinetic energy
        \alpha           | leaf area density [m]
        C_d              | tree canopy drag coefficient
        C_{p \epsilon 1} | coefficient with value 1.8 (Mochida et al. 2008)
    \endvartable

    Use with any variant of the k-\epsilon turbulence model.
    

    Reference:
        Mochida A., Tabata Y., Iwata T., Yoshino H., Examining tree canopy models for CFD prediction
        of wind environment at pedestrian level, J. Wind Eng. Ind.
        Aerodyn. 96 (2008) 1667-1677.
    
        Svensson, U., Häggkvist, K. A two-equation turbulence model for canopy flows. J. Wind Eng. Ind.
        Aerodyn. 35 (1990) 201–211.
    
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
