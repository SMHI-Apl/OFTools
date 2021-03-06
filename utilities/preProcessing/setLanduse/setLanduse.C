/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  Application
  setLanduse

  Description
  Set landuse values and heights

  \*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include <fstream>
#include "fileName.H"
#include "wallFvPatch.H"
#include "argList.H"
#include "groundDist.H"
#include "nutkAtmRoughWallFunctionFvPatchScalarField.H"
//Include classes for landuse
#include "Raster.H"

#ifndef namespaceFoam
#define namespaceFoam
using namespace Foam;
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool foundCode(List<tensor> codes, label code)
{
  forAll(codes,codei)
    {
      tensor codeTensor=codes[codei];
      if(codeTensor.xx()==code)
	return true;
    }
  return false;
}

label maxCode(const List<tensor> &codes)
{
  scalar maxCode=0;
  forAll(codes,codei)
    {
      tensor codeTensor=codes[codei];
      if(codeTensor.xx()>maxCode)
	maxCode=codeTensor.xx();
    }
  return label(maxCode);
}

void mapCodeIndices(const List<tensor> &codes,labelList &indexMap)
{
  forAll(codes,codei)
    {
    tensor codeTensor=codes[codei];
    indexMap[label(codeTensor.xx())]=codei;
    }
}

scalar integrateLAD(scalarList dist,scalar h,scalar maxLAD)
{
  //uses simpsons rule to integrate between the points
  scalar LAI=0;
  for(label fi=1;fi<dist.size();fi++)
    {
      if(dist[fi]==0 or dist[fi-1]==0)
	LAI+=1/2.0*abs(dist[1]-dist[0])*maxLAD*0.1*h;
      else
	LAI+=min(dist[fi],dist[fi-1])*maxLAD*0.1*h
	  +1/2.0*abs(dist[fi]-dist[fi-1])*0.1*h;
    }
  return LAI;
}

scalar getMaxLAD(scalarList dist,scalar h,scalar LAI)
{
  scalar maxLAD=0;
  scalar LAItmp=0;
  while(LAItmp<LAI)
    {
      LAItmp=integrateLAD(dist,h,maxLAD);
      maxLAD+=0.025;
    }
  return maxLAD;
}


// void setInlet(volVectorField &U_, scalar direction, scalar Umag) {
//   const scalar pi=3.1415926536;
//   vector Uref(0, 0 ,0);
//   scalar dirRad=pi/180*(90.-direction);
//   Uref.x()=-1*Foam::cos(dirRad);
//   Uref.y()=-1*Foam::sin(dirRad);
//   if(Uref.x()<1e-6 && Uref.x()>-1e-6)
//     Uref.x()=0.0;
//   if(Uref.y()<1e-6 && Uref.y()>-1e-6)
//     Uref.y()=0.0;

//   int number_of_patches=1;

//   const fvMesh & mesh=U_.mesh();
//   if (direction==360 || direction==0){
//     patchI = mesh.boundaryMesh().findPatchID("north");
//     Info<< "Found patch north at index " << patchI << endl;
//   }
//   else if (direction == 90) {
//     patchI = mesh.boundaryMesh().findPatchID("east");
//     Info<< "Found patch east at index " << patchI << endl;
//   }
//   else if (direction == 180){
//     patchI = mesh.boundaryMesh().findPatchID("south");
//     Info<< "Found patch south at index " << patchI << endl;
//   }
//   else if (direction == 270){
//     patchI = mesh.boundaryMesh().findPatchID("west");
//     Info<< "Found patch west at index " << patchI << endl;
//   }
//   else
//     {
//       number_of_patches=2;
//       if(direction>0 && direction <90)
// 	{
// 	  patchI=mesh.boundaryMesh().findPatchID("north");
// 	  patchII=mesh.boundaryMesh().findPatchID("east");
// 	  Info<< "Found patches north and east at index " << patchI <<" and" << patchII<< endl;
// 	}
//       else if(direction>90 && direction <180)
// 	{
// 	  patchI=mesh.boundaryMesh().findPatchID("east");
// 	  patchII=mesh.boundaryMesh().findPatchID("south");
// 	  Info<< "Found patches east and south at index " << patchI <<" and" << patchII<< endl;
// 	}
//       else if(direction>180 && direction <270)
// 	{
// 	  patchI=mesh.boundaryMesh().findPatchID("south");
// 	  patchII=mesh.boundaryMesh().findPatchID("west");
// 	  Info<< "Found patches south and west at index " << patchI <<" and" << patchII<< endl;
// 	}
//       else
// 	{
// 	  patchI=mesh.boundaryMesh().findPatchID("west");
// 	  patchII=mesh.boundaryMesh().findPatchID("north");
// 	  Info<< "Found patches west and north at index " << patchI <<" and" << patchII<< endl;
// 	}
//       if (patchII == -1)
// 	{
// 	  Info<< "Cannot find inlet patch\n" << endl;
// 	  FatalErrorIn(args.executable())
// 	    << "Cannot find inlet patch" << exit(FatalError);
// 	}
//     }

//   // Declaring and initializing references to boundary patches
//   const polyPatch& pp = mesh.boundaryMesh()[patchI];
//   fixedValueFvPatchVectorField& inletVelocity =
//     refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchI]);
//   fixedValueFvPatchScalarField& inletTurb =
//     refCast<fixedValueFvPatchScalarField>(k.boundaryField()[patchI]);
//   fixedValueFvPatchScalarField& inletEps =
//     refCast<fixedValueFvPatchScalarField>(epsilon.boundaryField()[patchI]);

//   Info<< "Wind mag. is " << Umag << " ,Wind dir. is "<< direction <<" degrees" << endl;
//   Info<< "X comp. of U at 10 m is " << U10.x() << endl;
//   Info<< "Y comp. of U at 10 m is " << U10.y() << endl;

// }



void setLanduse(label landuseCode, wordList sourcePatches,
		volScalarField &landuse_, volScalarField &LAD_,
                volScalarField &z0_, volScalarField &nut_,
                scalarList heightDist,
		List<tensor> &landuseList, Raster &lu,
		Switch readFromRaster)
{
  const fvMesh & mesh=landuse_.mesh();
  labelHashSet patchIDs(sourcePatches.size());

  // Create a mapping beween lu-code and index in tensorList
  // labelList indexMap(maxCode(landuseList));
  // mapCodeIndices(landuseList, indexMap);

  HashTable<label,label> indexMap(landuseList.size());
  forAll(landuseList,codei)
    {
      indexMap.insert(label(landuseList[codei].xx()),codei);
    }
  
  forAll(sourcePatches,sp)
    {
        //      Info<< "Processing patch1 " << sourcePatches[sp] << endl;
      label patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
      patchIDs.insert(patchI);
      const polyPatch& pp = mesh.boundaryMesh()[patchI];
      forAll(landuse_.boundaryField()[patchI],facei)
	{
          
	  scalar x=pp.faceCentres()[facei].x();
	  scalar y=pp.faceCentres()[facei].y();

	  if(readFromRaster)
              landuseCode=label(lu.getValue(double(x),double(y)));
              
          //Info<< "Set lu code " << scalar(landuseCode) << " in (x,y) = (" << x << ','  << y << ")" << endl;
	  landuse_.boundaryField()[patchI][facei]=scalar(landuseCode);
	}

      forAll(z0_.boundaryField()[patchI],facei)
        {
              
          //Info<< "Set lu code " << scalar(landuseCode) << " in (x,y) = (" << x << ','  << y << ")" << endl;
          scalar landuseCode = landuse_.boundaryField()[patchI][facei];
          tensor landuse_parameters = landuseList[indexMap[landuseCode]];
          z0_.boundaryField()[patchI][facei] = landuse_parameters.yy();; 
        }

      Foam::incompressible::nutkAtmRoughWallFunctionFvPatchScalarField& wallNut =
	refCast<Foam::incompressible::nutkAtmRoughWallFunctionFvPatchScalarField>(nut_.boundaryField()[patchI]);
      scalarField& z0 = wallNut.z0();

      forAll(landuse_.boundaryField()[patchI],facei)
	{
          scalar landuseCode = landuse_.boundaryField()[patchI][facei];
          tensor landuse_parameters = landuseList[indexMap[landuseCode]];
          z0[facei] =  landuse_parameters.yy();
          landuse_.boundaryField()[patchI][facei];
	}

    }

  volScalarField d
    (
     IOobject
     (
      "d",
     "constant",
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar("d", dimLength, SMALL),
     calculatedFvPatchScalarField::typeName
     );
  d = (const volScalarField&) groundDist(mesh,patchIDs).y();

  forAll(d.internalField(),celli)
    {
      scalar x=mesh.C()[celli].x();
      scalar y=mesh.C()[celli].y();
      // scalar z=mesh.C()[celli].z();

      if(readFromRaster)
	landuseCode=label(lu.getValue(double(x),double(y)));

      tensor code=landuseList[indexMap[landuseCode]];
      scalar height=code.yz();
      scalar groundDistance=d.internalField()[celli];
      if (groundDistance<height && height != 0)
	{
	  landuse_.internalField()[celli]=landuseCode;
	  label distIndex = label(round(groundDistance/height*10));//assumes that the heightDistribution is described by 11 values(step of 0.1*h) thus giving indices 0-10
	  scalar LADMax=code.zx();
	  LAD_.internalField()[celli]=heightDist[distIndex]*LADMax;
	}
    }
}


int main(int argc, char *argv[])
{
  // argList::validArgs.append(".asc file");
  // argList::validOptions.insert("raster", "rasterFile");
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readLanduseDict.H"

  label categoryNum = landuseList.size();

  if(categoryNum==0)
    FatalErrorIn(args.executable())
      << "Missing landuse classes in dictionary landuseData"
      << exit(FatalError);

  if(sourcePatches.size() == 0 && !readFromRaster)
    Info<< "Warning: No source patches given for landuse data";

  if((patchLanduse.size() != sourcePatches.size()) && !readFromRaster)
    FatalErrorIn(args.executable())
      << "Wrong number of rows in patchLanduse,"
      <<"should be equal to number of source patches "
      << exit(FatalError);

  Raster lu;
  if(readFromRaster)
    {
      if(!lu.read(rasterFileName.c_str()))
	FatalErrorIn(args.executable())
	  << "Cannot read file " << rasterFileName << exit(FatalError);

      lu.xll = lu.xll - double(subtractedX);
      lu.xur = lu.xur - double(subtractedX);
      lu.yll = lu.yll - double(subtractedY);
      lu.yur = lu.yur - double(subtractedY);

      Info << "-----------Raster specification-------------" << endl;
      Info << "Extent: " << lu.xll << "< X <" << lu.xur << " and "
	   << lu.yll << "< Y <" << lu.yur << endl;
      Info << "Dimensions: " << "cellsize=" << lu.cellsize
	   << " nrows= " << lu.nrows << " ncols= " << lu.ncols << endl;
    }

  //Check that all codes of raster exist in code table
  //   if(!foundCode(landuseList))
  // FatalErrorIn(args.executable()) <<": landuse code not found in landuseList in dictionary landuseData" <<exit(FatalError);

  //Initializing the landuse internalField and values on wall patches
  landuse.internalField() = -1;
  z0.internalField() = -1;
  const fvPatchList& patches = mesh.boundary();

  forAll(patches,patchI)
    {
      const fvPatch& curPatch = patches[patchI];
      if(isType<wallFvPatch>(curPatch))
	landuse.boundaryField()[patchI]=-1;
    }

    forAll(patches,patchI)
    {
      const fvPatch& curPatch = patches[patchI];
      if(isType<wallFvPatch>(curPatch))
	z0.boundaryField()[patchI]=-1;
    }

  forAll(landuseList,codei)
    {
      scalar LAI=landuseList[codei].xz();
      scalar height = landuseList[codei].yz();
      if(landuseList[codei].zx()==-1)
	{
	scalar maxLAD = getMaxLAD(heightDistribution,height,LAI);
	landuseList[codei].zx()=maxLAD;
	}
    }

  Info << "Landuse categories" << endl;
  Info << "id\tCd\tLAI\tfrac\tz0\theight\tmaxLAD" << endl;
  forAll(landuseList, codei)
    {
      tensor code = landuseList[codei];
      Info << code.xx() << "\t"<< code.xy() << "\t" << code.xz()
	   << "\t" << code.yx() << "\t" << code.yy() << "\t" << code.yz()
	   << "\t" << code.zx() << endl;
    }

  if (readFromRaster)
    {
      forAll(sourcePatches, sp)
	{
	  label patchI = mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
	  if (patchI == -1)
	    FatalErrorIn(args.executable())
	      << "Cannot find patch "<< sourcePatches[sp] << exit(FatalError);
	  else
	    Info<< "Found source patch " << sourcePatches[sp]
		<< " at patch with index " << patchI << endl;
	}

      //a dummy value is given for patchLanduseCode
      setLanduse(0, sourcePatches, landuse, LAD, z0,
                 nut, heightDistribution, landuseList, lu,
		 readFromRaster);
    }
  else
    {
      wordList currSourcePatches(1);
      forAll(sourcePatches,sp)
	{
	  label patchI = mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
	  if (patchI == -1)
	    FatalErrorIn(args.executable())
	      << "Cannot find patch " << sourcePatches[sp] << exit(FatalError);
	  else
	    Info<< "Found source patch " << sourcePatches[sp]
		<< " at patch with index " << patchI << endl;
	  currSourcePatches[0] = sourcePatches[sp];
	  Info << "Running setLanduse function" << endl;
	  setLanduse(patchLanduse[sp], currSourcePatches,
            landuse,LAD, z0, nut, heightDistribution,
            landuseList, lu,readFromRaster);
	}
    }

    landuse.write();
    LAD.write();
    nut.write();
    z0.write();
  Info<< "setLanduse finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //

