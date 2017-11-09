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

#include "DalpeMassonCanopySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"
#include "groundDist.H"
#include "Raster.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(DalpeMassonCanopySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        DalpeMassonCanopySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::DalpeMassonCanopySource::DalpeMassonCanopySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
  option(name, modelType, dict, mesh),
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


void Foam::fv::DalpeMassonCanopySource::checkData() const
{
  if (betaP_.value() <= 0)
    {
        FatalErrorInFunction
           << "betaP must be greater than zero"
           << exit(FatalIOError);
    }
  if (betaD_.value() <= 0)
    {
        FatalErrorInFunction
           << "betaD must be greater than zero"
           << exit(FatalIOError);
    }
  if (C4_.value() <= 0)
    {
        FatalErrorInFunction
           << "C4 must be greater than zero"
           << exit(FatalIOError);
    }
  if (C5_.value() <= 0)
    {
        FatalErrorInFunction
           << "C5 must be greater than zero"
           << exit(FatalIOError);
    }
  
  if(landuseTable_.size() == 0)
    FatalErrorInFunction
      << "Missing landuse classes for canopySource in constant/fvOptions"
      << exit(FatalError);

  if(sourcePatches_.size() == 0 && !readFromRaster_)
    FatalErrorInFunction
      << "No source patches given for landuse data"
      << exit(FatalError);
}


void Foam::fv::DalpeMassonCanopySource::addSup
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


void Foam::fv::DalpeMassonCanopySource::addSup
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


void Foam::fv::DalpeMassonCanopySource::addSup
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
     betaP_*canopy*pow(mag(U),3) - fvm::SuSp(betaD_*canopy*mag(U), k)
    );

    eqn +=  Sk;
  }
  else if (eqn.psi().name() == word("epsilon")) {

    const volScalarField& epsilon = eqn.psi();
    const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    
    fvMatrix<scalar> Sepsilon
    (
     fvm::SuSp(canopy/k*(C4_*betaP_*pow(mag(U),3) - C5_*betaD_*k*mag(U)), epsilon)
    );

    eqn += Sepsilon;    
  }
}


void Foam::fv::DalpeMassonCanopySource::addSup
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


bool Foam::fv::DalpeMassonCanopySource::read(const dictionary& dict)
{
    if (!option::read(dict))
      return false;

    
      if (active()) {
        coeffs_.readIfPresent("betaP", betaP_);
        coeffs_.readIfPresent("betaD", betaD_);
        coeffs_.readIfPresent("C4", C4_);
        coeffs_.readIfPresent("C5", C5_);
        coeffs_.lookup("sourcePatches") >> sourcePatches_;
        rasterOriginX_ = coeffs_.lookupOrDefault("rasterOriginX", 0, false, false);
        rasterOriginY_ = coeffs_.lookupOrDefault("rasterOriginY", 0, false, false);
        readFromRaster_ = coeffs_.lookupOrDefault("readFromRaster", false, false, false);
        writeFields_ = coeffs_.lookupOrDefault("writeFields", false, false, false);
        
        // read landuse raster
        if (readFromRaster_ ) {
          coeffs_.lookup("rasterFileName") >>  rasterFileName_;
          
          if(!raster_.read(rasterFileName_.c_str()))
            FatalErrorInFunction
              << "Cannot read file " << rasterFileName_
              << exit(FatalError);

          raster_.xll = raster_.xll - double(rasterOriginX_);
          raster_.xur = raster_.xur - double(rasterOriginX_);
          raster_.yll = raster_.yll - double(rasterOriginY_);
          raster_.yur = raster_.yur - double(rasterOriginY_);
          
          Info << "    -----------------------Raster specification----------" << endl
               << "    Extent: " << raster_.xll << "< X <" << raster_.xur << " and "
               << "    " << raster_.yll << "< Y <" << raster_.yur << endl
               << "    Dimensions: " << " nrows= " << raster_.nrows << " ncols= " << raster_.ncols << endl
               << "    Cellsize=" << raster_.cellsize << endl;
        }
          
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

        if (!readFromRaster_) {
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
              FatalIOErrorInFunction(
                "void Foam::fv::DalpeMassonCanopySource::read("
                "const dictionary&)"
              ) << "Cannot find patch " << sourcePatches_[i] << exit(FatalIOError);

            label patchLanduseCode = patchLanduseList[i];
            landuseClass lu = landuseTable_[patchLanduseCode];
            
            patchLanduseTable_.insert(patchID, landuseTable_[lu.code()]);
          }
        }
        checkData();
        calculateCanopy();
      }
      return true;
}

void Foam::fv::DalpeMassonCanopySource::setPatchLanduse(
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
      if (readFromRaster_) {
        scalar x=pp.faceCentres()[facei].x();
        scalar y=pp.faceCentres()[facei].y();
        lu = landuseTable_[label(raster_.getValue(double(x),double(y)))];
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

  forAll(d.internalField(),celli)
    {

      // get landuse class
      landuseClass lu;
      if (readFromRaster_) {

        scalar x=mesh.C()[celli].x();
        scalar y=mesh.C()[celli].y();
        lu = landuseTable_[label(raster_.getValue(double(x),double(y)))];
      }
      else {
        lu = patchLanduseTable_[patch];
      }

      // set landuse tp canopy height and LAD according to LADProfile
      scalar patchDistance = d.internalField()[celli];
      if (patchDistance < lu.height() && lu.height() != 0)
	{
	  landuse.primitiveFieldRef()[celli]=scalar(lu.code());
          scalar deltaz = lu.height() / scalar(lu.LADProfile().size());
          label LADProfileValueIndex = label(round(patchDistance / deltaz));
	  LAD.primitiveFieldRef()[celli] = lu.LADProfile()[LADProfileValueIndex] * lu.LADmax();
	}
    }
}


void Foam::fv::DalpeMassonCanopySource::calculatePatchDistance(label patch, volScalarField& d)
{

  labelHashSet sourcePatchIDs(1);  
  sourcePatchIDs.insert(patch);
  // calculate distance to sourcePatches
  d = (const volScalarField&) groundDist(mesh_,sourcePatchIDs).y();
  return;
}

void Foam::fv::DalpeMassonCanopySource::calculateCanopy()
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
         "void Foam::fv::DalpeMassonCanopySource::setPatchLanduse("
         "landuseClass,"
         "label,"
         "volScalarField&,"
         "volScalarField&,"
         "volScalarField&)",
         coeffs_
      ) << "Cannot find patch " << patchName << exit(FatalIOError);

    if (mesh_.boundaryMesh()[patchName].size() == 0)
      FatalIOErrorIn(
         "void Foam::fv::DalpeMassonCanopySource::setPatchLanduse("
         "landuseClass,"
         "label,"
         "volScalarField&,"
         "volScalarField&,"
         "volScalarField&)",
         coeffs_
      ) << "Source patch " << patchName << " has zero faces"<< exit(FatalIOError);
    
    calculatePatchDistance(patchID, d);
    d.write();
    
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
