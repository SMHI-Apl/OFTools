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
    David Segersson, Swedish Meteorological and Hydrological Institute
    david.segersson@smhi.se

Description
    Base class for momentum and turbulence source/sink-terms for tree canopy.


Usage
    Example usage:
    \verbatim
    canopy
    {
        type                       canopySource;
        active                     on;              // toggle the option
        writeFields                on;              // write to disk
        readLanduseFromRaster      on;              // read landuse from raster
        readCanopyHeightFromRaster on;              // read canopy height from raster
        sourcePatches              (ground forest); // for which patches to set landuse
        patchLanduse               (0 1);           // landuse code per patch (if not read from raster)

        landuse
        {

            low_birch                               // name of landuse class
            {
                code 1;                             // code of landuse class
                Cd 0.2;                             // drag coefficient of trees
                LAI 2.15;                           // Leaf Area Index (used if LADmax is not specified)
                LADmax 1.2;                         // Maximum Leaf Area Density [m ^-1]
                z0 0.06;                            // roughness length [m]
                height 7.5;                         // tree canopy height (if not read from raster)
                
                // Vertical profile of Leaf Area Density
                // first value is closest to ground
                // each value represents an equal share of the tree height
                // e.g. for a 4 m tree and 4 values, each value will represent 1 m
                // values are scaled so that the highest will correspond to LADMax                

                LADProfile ( 0.05 0.1 0.15 0.35 1.1 0.9 0.5 0.2 0.15 0.05 0.01 );
            };

            grass
            {
                code 0;
                Cd 0.2;
                LAI 0;
                z0 0.06;
                height 0;
                LADProfile ( 0.05 0.1 0.15 0.35 1.1 0.9 0.5 0.2 0.15 0.05 0.01 );
            };

        };
    }
    \endverbatim
    
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
        autoPtr<volScalarField> canopy_;

        HashTable<landuseClass, label> landuseTable_;
        HashTable<landuseClass, label> patchLanduseTable_;

        // read landuse from raster instead of specifying per patch
        Switch readLanduseFromRaster_;
        Raster landuseRaster_;

        // read canopy height from raster instead of specifying per landuse class
        Switch readCanopyHeightFromRaster_;
        Raster canopyHeightRaster_;

        // translation from coordiantes in raster to coordinates in mesh
        scalar originOffsetX_;
        scalar originOffsetY_;

        Switch writeFields_;

        // if variables are present on disk, they are read instead of calculated
        // this allows setting the data with external programs
        bool _LAD_from_disk;
        bool _landuse_from_disk;
        bool _z0_from_disk;

        //- Source terms to momentum equation (for solvers with and without explicit density)
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

        // Member Functions
        Raster readRaster(fileName rasterPath);
        void calculatePatchDistance(label patch, volScalarField &d);
        void setPatchLanduse(label patch, volScalarField &landuse, 
                             volScalarField &LAD, volScalarField &z0,
                             volScalarField &nut, volScalarField &d);
        void calculateCanopy();
  void readLanduseClasses();

  
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
