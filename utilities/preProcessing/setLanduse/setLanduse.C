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

void setLanduse(label landuseCode, wordList sourcePatches, volScalarField &landuse_, volScalarField &LAD_, volScalarField &nut_, scalarList heightDist, List<tensor> &landuseList, Raster &lu, word dataSource)
{
  const fvMesh & mesh=landuse_.mesh();
  labelHashSet patchIDs(sourcePatches.size());
  forAll(sourcePatches,sp)
    {
      label patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
      patchIDs.insert(patchI);

      const polyPatch& pp = mesh.boundaryMesh()[patchI];
      //      fixedValueFvPatchScalarField& patchLanduse=
      //	refCast<fixedValueFvPatchScalarField>(landuse_.boundaryField()[patchI]);
      forAll(landuse_.boundaryField()[patchI],facei)
	{
	  scalar x=pp.faceCentres()[facei].x();
	  scalar y=pp.faceCentres()[facei].y();
	  if(dataSource=="fromFile")
	    landuseCode=label(lu.getValue(double(x),double(y)));
	  //patchLanduse[facei]=landuseCode;
	  landuse_.boundaryField()[patchI][facei]=scalar(landuseCode);
	}

      Foam::incompressible::nutkAtmRoughWallFunctionFvPatchScalarField& wallNut =
	refCast<Foam::incompressible::nutkAtmRoughWallFunctionFvPatchScalarField>(nut_.boundaryField()[patchI]);
      scalarField& z0 = wallNut.z0();
      forAll(landuse_.boundaryField()[patchI],facei)
	{
	  z0[facei] = landuse_.boundaryField()[patchI][facei];
	}

    }
  // Create a mapping beween lu-code and index in tensorList
  // labelList indexMap(maxCode(landuseList));
  // mapCodeIndices(landuseList, indexMap);

  HashTable<label,label> indexMap(landuseList.size());
  forAll(landuseList,codei)
    {
      indexMap.insert(label(landuseList[codei].xx()),codei);
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
      scalar z=mesh.C()[celli].z();
      
      if(dataSource=="fromFile")
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
  argList::validOptions.insert("raster", "rasterFile");  
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readLanduseDict.H"


  fileName rasterName;
  bool readRaster = args.options().found("raster");
  word dataSource("Manual");
  if (readRaster)
    {
      rasterName=fileName(args.options()["raster"]);
      dataSource="fromFile";
    }

  label categoryNum = landuseList.size();

  if(categoryNum==0)
    FatalErrorIn(args.executable()) 
      << "Missing landuse classes in dictionary landuseData"
      << exit(FatalError);

  if(sourcePatches.size()==0 && dataSource=="Manual")
    Info<< "Warning: No source patches given for landuse data";

  if(patchLanduse.size()!=sourcePatches.size() && dataSource == "Manual")
    FatalErrorIn(args.executable())
      << "Wrong number of rows in patchLanduse,"
      <<"should be equal to number of source patches "
      << exit(FatalError);
  
  Raster lu;
  if(dataSource=="fromFile")
    {
      if(!lu.read(rasterName.c_str()))
	FatalErrorIn(args.executable())
	  << "Cannot read file " << rasterName << exit(FatalError);
      
      lu.xll=lu.xll-double(subtractedX);
      lu.xur=lu.xur-double(subtractedX);
      lu.yll=lu.yll-double(subtractedY);
      lu.yur=lu.yur-double(subtractedY);
        
      Info<<"-----------Raster specification-------------"<<endl;
      Info<<"Extent: "<<lu.xll<<"< X <"<<lu.xur<<" and "
	  <<lu.yll<<"< Y <"<<lu.yur<<endl;
      Info<<"Dimensions: "<<"cellsize="<<lu.cellsize
	  <<" nrows= "<<lu.nrows<<" ncols= "<<lu.ncols<<endl; 
    }      
 
  //Check that all codes of raster exist in code table
  //   if(!foundCode(landuseList))
  // FatalErrorIn(args.executable()) <<": landuse code not found in landuseList in dictionary landuseData" <<exit(FatalError);
   
  //Initializing the landuse internalField and values on wall patches
  landuse.internalField()=-1;
  const fvPatchList& patches=mesh.boundary();
  forAll(patches,patchI)
    {
      const fvPatch& curPatch=patches[patchI];
      if(isType<wallFvPatch>(curPatch))
	landuse.boundaryField()[patchI]=-1;
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

  Info<<"Landuse categories"<<endl;
  Info<<"id\tCd\tLAI\tfrac\tz0\theight\tmaxLAD"<<endl;
  forAll(landuseList,codei)
    {
      tensor code=landuseList[codei];
      Info<<code.xx()<<"\t"<<code.xy()<<"\t"<<code.xz()
	  <<"\t"<<code.yx()<<"\t"<<code.yy()<<"\t"<<code.yz()
	  <<"\t"<<code.zx()<<endl;
    }

  if (dataSource=="fromFile")
    {
      forAll(sourcePatches,sp)
	{
	  label patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
	  if (patchI == -1)
	    FatalErrorIn(args.executable())
	      << "Cannot find patch "<< sourcePatches[sp] << exit(FatalError);
	  else
	    Info<< "Found source patch "<< sourcePatches[sp]
		<<" at patch with index " << patchI << endl;
	}
      
      setLanduse(0,sourcePatches,landuse,LAD, nut, heightDistribution,landuseList,lu,dataSource);//a dummy value is given for patchLanduseCode 
    }
  else
    {
      wordList currSourcePatches(1);
      forAll(sourcePatches,sp)
	{
	  label patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
	  if (patchI == -1)
	    FatalErrorIn(args.executable())
	      << "Cannot find patch "<< sourcePatches[sp] << exit(FatalError);
	  else
	    Info<< "Found source patch "<< sourcePatches[sp]
		<<" at patch with index " << patchI << endl;
	  currSourcePatches[0]=sourcePatches[sp];
	  Info<<"Running setLanduse function"<<endl;
	  setLanduse(patchLanduse[sp],currSourcePatches,landuse,LAD, nut, heightDistribution,landuseList,lu,dataSource);
	}
    }
  
    landuse.write();
    LAD.write();
    nut.write();
  Info<< "setLanduse finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //

