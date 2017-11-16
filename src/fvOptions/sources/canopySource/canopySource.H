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
    Foam::fv::canopySource

Group
    grpFvOptionsSources

Author
    David Segersson, Swedish Meteorological and Hydrological Institute, david.segersson@smhi.se

Description
    Base class for momentum and turbulence source/sink-terms for tree canopy.
    
SourceFiles
    canopySource.C

\*---------------------------------------------------------------------------*/

#ifndef canopySource_H
#define canopySource_H

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
                 Class canopySource Declaration
\*---------------------------------------------------------------------------*/

class canopySource
:
    public option
{

protected:

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

        //- Source term to momentum equation
        void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Disallow default bitwise copy construct
        canopySource(const canopySource&);

        //- Disallow default bitwise assignment
        void operator=(const canopySource&);

        void checkData() const;

        // Private Member Functions
        void calculatePatchDistance(label patch, volScalarField &d);
        void setPatchLanduse(label patch, volScalarField &landuse, 
                             volScalarField &LAD, volScalarField &z0,
                             volScalarField &nut, volScalarField &d);
        void calculateCanopy();

  
public:


    //- Runtime type information
    TypeName("canopySource");


    // Constructors

        //- Construct from components
        canopySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~canopySource()
    {}


    // Member Functions

        //- Read dictionary
        bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
