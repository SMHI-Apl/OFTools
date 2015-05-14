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
  createArenaFaceSets

  Description
  Creates face sets corresponding to createFaceSetDict

  \*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include <fstream>
#include "fileName.H"
#include "wallFvPatch.H"
#include "faceSet.H"
#include "PtrList.H"
#include <iostream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readFaceSetsDict.H"

  // Finding indexes to given patches 

  Info<<"Creating face sets"<<endl;
  PtrList<faceSet> faceSets((xlims.size()-1)*(ylims.size()-1)*(zlims.size()-1));
  List<scalar> numFaces((xlims.size()-1)*(ylims.size()-1)*(zlims.size()-1));
  
  for(label i=0;i<numFaces.size();i++)
    {
      numFaces[i]=0;
    }
  

  // Declaring and initializing references to boundary patches

  label maxPatchSize=0;
  forAll(patchNames,namei)
    {
      label patchI;
      patchI=mesh.boundaryMesh().findPatchID(patchNames[namei]);

      if (patchI == -1)
	{
	  FatalErrorIn(args.executable())
	    << "Cannot find patch: "<<patchNames[namei]<<exit(FatalError);
	}

      label patchSize=mesh.boundaryMesh()[patchI].faceCentres().size();
  
      if (patchSize > maxPatchSize)
	maxPatchSize=patchSize;
    }

  for(label i=0;i<faceSets.size();i++)
    {
      std::ostringstream outStr;
      outStr<<prefix<<"_"<<i+1;
      word setName(outStr.str());
      faceSets.set(i,new faceSet(mesh,setName,maxPatchSize));
    }	
  
  forAll(patchNames,namei)
    {

      label patchI;
      patchI=mesh.boundaryMesh().findPatchID(patchNames[namei]);
  
      if (patchI == -1)
	{
	  FatalErrorIn(args.executable())
	    << "Cannot find patch: "<<patchNames[namei]<<exit(FatalError);
	}
      const polyPatch& pp = mesh.boundaryMesh()[patchI];
  
      Info<<"Looping over faces of patch: "<<patchNames[namei]<<endl;      // Looping over all faces of the chosen patch
      forAll(pp.faceCentres(), faceI)
	{  
	  //get coordinates of face centre
	  const point& faceCentre= pp.faceCentres()[faceI];
	
	  scalar x=faceCentre.x();
	  scalar y=faceCentre.y();
	  scalar z=faceCentre.z();
	  vector faceNormal=pp.faceAreas()[faceI];
	  scalar theta=0;
          scalar pi=3.1415926535897932;

          // calculate theta, angle between normal and faceNormal
	  theta = Foam::acos( (normal & faceNormal)/(mag(normal)*mag(faceNormal)) );
           
          // converting from rad to grad
	  theta=theta*180/pi;
         
	  label faceIndex=pp.start()+faceI;
	  label setIndex=0;
          label xyIndex=0;

	  if(theta <= tol)
	    {

	      label yi=1;
	      bool foundY=false;
	      for(yi=1;yi<ylims.size();yi++)
		{
		  scalar ymax=ylims[yi];
		  if(y<ymax)
		    {
		      foundY=true;
		      break;
		    }
		}
	      label xi=1;
	      bool foundX=false;
	      for(xi=1;xi<xlims.size();xi++)
		{
		  scalar xmax=xlims[xi];
		  if(x<xmax)
		    {
		      foundX=true;
		      break;
		    }
		}		      
              label zi=1;
	      bool foundZ=false;
	      for(zi=1;zi<zlims.size();zi++)
		{
		  scalar zmax=zlims[zi];
		  if(z<zmax)
		    {
		      foundZ=true;
		      break;
		    }
		}		      
	      if( !(foundX and foundY and foundZ))
		Info<<"Outer limits for face sets are smaller than patch extent"<<FatalErrorIn(args.executable());
	  
	      xyIndex=(yi-1)*(xlims.size()-1)+xi-1;
              setIndex=(zi-1)*(xlims.size()-1)*(ylims.size()-1)+xyIndex;

	      faceSets[setIndex].insert(faceIndex);
	      numFaces[setIndex]=numFaces[setIndex]+1;
	    }
        }
    }
      for(label i=0;i<faceSets.size();i++)
	{
	  faceSets[i].write();
	}
      Info<<" Number of cell faces in faceSets:"<<endl;
      for(label seti=0;seti<numFaces.size();seti++)
	{
	  Info<<faceSets[seti].name()<<"\t"<<numFaces[seti]<<endl;
	}
      
      Info<< "createFaceSets finished successfully\n" << endl;
      return(0);
   
}
// ************************************************************************* //

