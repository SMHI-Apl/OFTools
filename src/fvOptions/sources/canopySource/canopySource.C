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

#include "canopySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"
#include "groundDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(canopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        canopySource,
        dictionary
    );

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::canopySource::canopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
  option(name, modelType, dict, mesh)
{
 
  if (active()) {
    
    fieldNames_.setSize(3);
    fieldNames_[0] = word("U");
    fieldNames_[1] = word("k");
    fieldNames_[2] = word("epsilon");
    applied_.setSize(fieldNames_.size(), false);

    read(dict);
  }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::canopySource::addSup
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


void Foam::fv::canopySource::addSup
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


void Foam::fv::canopySource::checkData() const
{

  if(landuseTable_.size() == 0)
    FatalErrorInFunction
      << "Missing landuse classes for canopySource in constant/fvOptions"
      << exit(FatalError);

  if(sourcePatches_.size() == 0 && !readLanduseFromRaster_)
    FatalErrorInFunction
      << "No source patches given for landuse data"
      << exit(FatalError);
}

Raster Foam::fv::canopySource::readRaster(fileName rasterPath){

        originOffsetX_ = coeffs_.lookupOrDefault("rasterOriginX", 0, false, false);
        originOffsetY_ = coeffs_.lookupOrDefault("rasterOriginY", 0, false, false);

        Raster raster;
        // read raster
        if(!raster.read(rasterPath.c_str()))
          FatalErrorInFunction
            << "Cannot read file " << rasterPath
            << exit(FatalError);
        raster.translate(double(originOffsetX_), double(originOffsetY_));
        
        Info << "    Raster: " << rasterPath << endl
             << "        Extent: " << raster.xll << "< X <" << raster.xur << " and "
             << "        " << raster.yll << "< Y <" << raster.yur << endl
             << "        Dimensions: " << " nrows= " << raster.nrows << " ncols= " << raster.ncols << endl
             << "        Cellsize=" << raster.cellsize << endl;
        return raster;
}

void Foam::fv::canopySource::readLanduseClasses()
{
  // get landuse definitions
  dictionary landuseDict = coeffs_.subDict("landuse");    
  wordList landuseNames(landuseDict.toc());
  
  Info << "    ---------------------Landuse categories--------------------" << endl;
  Info << "    name\tcode\tCd\tLAI\tz0\theight\tLADmax" << endl;
  forAll(landuseNames, i) {
      word name = landuseNames[i];
      //dictionary luDict = landuseDict.subDict(name);
      landuseClass lu(landuseDict, name);
      landuseTable_.insert(lu.code(), lu);
      Info << "    " << lu.name() << "\t"<< lu.code() << "\t" << lu.Cd()
           << "\t" << lu.LAI() << "\t" << lu.z0() << "\t" << lu.height()
           << "\t" << lu.LADmax() << endl;
  }

}


bool Foam::fv::canopySource::read(const dictionary& dict)
{
    if (!option::read(dict))
      return false;

    // if option is not active, it will not be read
    if (! active())
      return true;

    // for which patches to set landuse
    coeffs_.lookup("sourcePatches") >> sourcePatches_;

    readLanduseFromRaster_ = coeffs_.lookupOrDefault("readLanduseFromRaster", false, false, false);
    readCanopyHeightFromRaster_ = coeffs_.lookupOrDefault("readCanopyHeightFromRaster", false, false, false);
    writeFields_ = coeffs_.lookupOrDefault("writeFields", false, false, false);
    
    if (readLanduseFromRaster_ ) {
      fileName landuseRasterFileName;
      coeffs_.lookup("landuseRasterFileName") >>  landuseRasterFileName;
      landuseRaster_ = readRaster(landuseRasterFileName);
    }

    if (readCanopyHeightFromRaster_ ) {
      fileName canopyHeightRasterFileName;
      coeffs_.lookup("canopyHeightRasterFileName") >>  canopyHeightRasterFileName;
      canopyHeightRaster_ = readRaster(canopyHeightRasterFileName);
    }

    readLanduseClasses();

    // read landuse codes for each patch
    if (!readLanduseFromRaster_) {
      labelList patchLanduseList;
      coeffs_.lookup("patchLanduse") >> patchLanduseList;
      if((patchLanduseList.size() != sourcePatches_.size()))
        FatalErrorInFunction
          << "Wrong number of rows in patchLanduse,"
          <<"should be equal to number of source patches "
          << exit(FatalError);
          
      forAll(sourcePatches_, i) {
        label patchID = mesh_.boundaryMesh().findPatchID(sourcePatches_[i]);
        if (patchID == -1)
          FatalErrorInFunction
            << "Cannot find landuse source patch: "
            << sourcePatches_[i]
            << exit(FatalError);

        label patchLanduseCode = patchLanduseList[i];
        landuseClass lu = landuseTable_[patchLanduseCode];        
        patchLanduseTable_.insert(patchID, landuseTable_[lu.code()]);
      }
    }
    checkData();
    calculateCanopy();
    return true;
}

void Foam::fv::canopySource::setPatchLanduse(
                label patch,
                volScalarField &landuse,
                volScalarField &LAD,
                volScalarField &z0,
                volScalarField &nut,
                volScalarField &d)
{

  const fvMesh & mesh=landuse.mesh();
  const polyPatch& pp = mesh.boundaryMesh()[patch];


  Foam::nutkAtmRoughWallFunctionFvPatchScalarField& wallNut =
    refCast<Foam::nutkAtmRoughWallFunctionFvPatchScalarField>(nut.boundaryFieldRef()[patch]);

  scalarField& nutZ0 = wallNut.z0();
   
  Foam::fvPatchScalarField& patchLanduse = landuse.boundaryFieldRef()[patch];
  Foam::fvPatchScalarField& patchZ0 = z0.boundaryFieldRef()[patch];

  forAll(patchLanduse,facei)
    {
      
      landuseClass lu;
      if (readLanduseFromRaster_) {
        scalar x=pp.faceCentres()[facei].x();
        scalar y=pp.faceCentres()[facei].y();
        lu = landuseTable_[label(landuseRaster_.getValue(double(x),double(y)))];
      }
      else {
        lu = patchLanduseTable_[patch];
      }

      if (!_landuse_from_disk)
        patchLanduse[facei] = scalar(lu.code());
      
      if (!_z0_from_disk)
        patchZ0[facei] = lu.z0();

      // z0 for actual calculations is stored as scalarField 
      // within boundary condition dictionary of nut
      nutZ0[facei] = patchZ0[facei];
    }

  forAll(d.internalField(), celli)
    {

      // get landuse class
      landuseClass lu;

      scalar x = mesh.C()[celli].x();
      scalar y = mesh.C()[celli].y();
      scalar height = 0;

      if (readLanduseFromRaster_)
        lu = landuseTable_[label(landuseRaster_.getValue(double(x),double(y)))];
      else
        lu = patchLanduseTable_[patch];

      if (readCanopyHeightFromRaster_)
        height = canopyHeightRaster_.getValue(double(x), double(y));
      else
        height = lu.height();

      // set landuse up to canopy height and LAD according to LADProfile
      scalar patchDistance = d.internalField()[celli];
      if (patchDistance < height && height > 0)
	{
	  landuse.primitiveFieldRef()[celli] = scalar(lu.code());
	  LAD.primitiveFieldRef()[celli] = lu.LAD(patchDistance, height);
	}
    }
}


void Foam::fv::canopySource::calculatePatchDistance(label patch, volScalarField& d)
{

  labelHashSet sourcePatchIDs(1);  
  sourcePatchIDs.insert(patch);
  // calculate distance to sourcePatches
  d = (const volScalarField&) groundDist(mesh_,sourcePatchIDs).y();
  return;
}

void Foam::fv::canopySource::calculateCanopy()
{


  autoPtr<volScalarField> landuse;
  autoPtr<volScalarField> LAD;
  autoPtr<volScalarField> z0;

  IOobject landuseHeader
  (
      "landuse",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
  );

  if (landuseHeader.typeHeaderOk<volScalarField>(true)) {
    _landuse_from_disk = true;
    landuse.reset
    (
        new volScalarField
        (
            IOobject
            (
                "landuse",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
     );
  }
  else {
    landuse.reset
    (
        new volScalarField
        (
             IOobject
             (
                 "landuse",
                 mesh_.time().timeName(),
                 mesh_,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             mesh_,
             dimensionedScalar("landuse", dimensionSet(0,0,0,0,0,0,0), -1),
             calculatedFvPatchScalarField::typeName
        )
    );
  }

  IOobject z0Header
  (
      "z0",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
  );

  if (z0Header.typeHeaderOk<volScalarField>(true)) {
    _z0_from_disk = true;
    z0.reset
    (
        new volScalarField
        (
            IOobject
            (
                "z0",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
     );
  }
  else {
    z0.reset
    (
        new volScalarField
        (
             IOobject
             (
                 "z0",
                 mesh_.time().timeName(),
                 mesh_,
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             mesh_,
             dimensionedScalar("z0", dimLength, -1),
             calculatedFvPatchScalarField::typeName
        )
    );
  }

  IOobject LADHeader
  (
      "LAD",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
  );

  if (LADHeader.typeHeaderOk<volScalarField>(true)) {
    _LAD_from_disk = true;
    
    LAD.reset
    (
       new volScalarField
       (
           IOobject
           (
               "LAD",
               mesh_.time().timeName(),
               mesh_,
               IOobject::MUST_READ,
               IOobject::AUTO_WRITE
           ),
           mesh_
       )
    );
  }
  else 
  {
    LAD.reset
    (
        new volScalarField
        (
           IOobject
           (
               "LAD",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           mesh_,
           dimensionedScalar("lad", dimensionSet(0,-1,0,0,0,0,0), 0),
           calculatedFvPatchScalarField::typeName
        )
    );
  }

  // patch distance
  volScalarField d
  (
     IOobject
     (
          "d",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
     ),
     mesh_,
     dimensionedScalar("d",dimLength, SMALL),
     calculatedFvPatchScalarField::typeName
  );
  
  canopy_.reset
  (
     new volScalarField(
         IOobject
         (
            "canopy",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         mesh_,
         dimensionedScalar("canopy",dimensionSet(0,-1,0,0,0,0,0),0)
     )
  );

  volScalarField& nut = mesh_.lookupObjectRef<volScalarField>("nut");
  scalarField& canopy = canopy_->primitiveFieldRef();

  // set landuse, LAD and z0
  Info<<"\n    Setting landuse for patch: " << endl;
  forAll(sourcePatches_, sp) {
    word patchName = sourcePatches_[sp];
    Info<<"    -- " << patchName << endl; 
    label patchID = mesh_.boundaryMesh().findPatchID(patchName);
    
    if (patchID == -1)
      FatalIOErrorIn(
         "void Foam::fv::canopySource::setPatchLanduse("
         "landuseClass,"
         "label,"
         "volScalarField&,"
         "volScalarField&,"
         "volScalarField&)",
         coeffs_
      ) << "Cannot find patch " << patchName << exit(FatalIOError);


      calculatePatchDistance(patchID, d);    
      setPatchLanduse(patchID, landuse(), LAD(), z0(), nut, d);
  }

  // calculate canopy (Cd * LAD)
  forAll(landuse->internalField(),celli) {
    label landuseCode(landuse->internalField()[celli]);
    if (landuseCode != -1) {
      landuseClass lu = landuseTable_[landuseCode];
      canopy[celli] = lu.Cd() * LAD->internalField()[celli];
    }
  }


  if (writeFields_) {
    landuse->write();
    LAD->write();
    z0->write();
  }     
}


// ************************************************************************* //
