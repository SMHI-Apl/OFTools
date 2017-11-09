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
    Foam::fv::DalpeMassonCanopySource

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
        type            DalpeMassonCanopySource;
        active          on;
    }
    \endverbatim

SourceFiles
    DalpeMassonCanopySource.C

\*---------------------------------------------------------------------------*/

#ifndef DalpeMassonCanopySource_H
#define DalpeMassonCanopySource_H

#include "fvOption.H"
#include "dimensionedTypes.H"
#include "DynamicList.H"
#include "landuseClass.H"
#include "Raster.H"
#include "nutkAtmRoughWallFunctionFvPatchScalarField.H"
#include "groundDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class DalpeMassonCanopySource Declaration
\*---------------------------------------------------------------------------*/

class DalpeMassonCanopySource
:
    public option
{

private:

        // Canopy coefficients
        dimensionedScalar betaP_;
        dimensionedScalar betaD_;
        dimensionedScalar C4_;
        dimensionedScalar C5_;
  
        // member data for landuse
        wordList sourcePatches_;
        fileName rasterFileName_;
        scalar rasterOriginX_;
        scalar rasterOriginY_;
        Switch readFromRaster_;
        autoPtr<volScalarField> canopy_;
        HashTable<landuseClass, label> landuseTable_;
        HashTable<landuseClass, label> patchLanduseTable_;
        Raster raster_;
        Switch writeFields_;
        bool _LAD_from_disk;
        bool _landuse_from_disk;
        bool _z0_from_disk;

        //- Disallow default bitwise copy construct
        DalpeMassonCanopySource(const DalpeMassonCanopySource&);

        //- Disallow default bitwise assignment
        void operator=(const DalpeMassonCanopySource&);

        void checkData() const;

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
        void calculatePatchDistance(label patch, volScalarField &d);
        void setPatchLanduse(label patch, volScalarField &landuse, 
                             volScalarField &LAD, volScalarField &z0,
                             volScalarField &nut, volScalarField &d);
        void calculateCanopy();

  
public:


    //- Runtime type information
    TypeName("DalpeMassonCanopySource");


    // Constructors

        //- Construct from components
        DalpeMassonCanopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~DalpeMassonCanopySource()
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
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to scalar equations (k & epsilon)
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Source term to compressible scalar equation (k & epsilon)
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read dictionary
        virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
