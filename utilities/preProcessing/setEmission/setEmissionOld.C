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
  setEmission

  Description
  Set emissions from traffic in a volume up to specified height above the road.
  The emission is normalized so that it corresponds to a emission of 1 g/s for each emissionGroup.
  The separate roads (patches) within the emissionGroup are weighted
  To get other emissions, multiply the concentration field for each emissionGroup with the real emission in g/s

  \*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "fileName.H"
#include "speciesTable.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readEmissions.H" 
#   include "createFields.H"
  // * ** * * * * * * * * * * * * * * * * * * * * * * //  
  
  //Checking for errors in input
  //------------------------------------------------------------------
  Info<<"Started setEmission"<<endl;

  if(emissionGroups.size()==0)
    FatalErrorIn(args.executable()) << "Road source patches not specified"<< exit(FatalError);

 
  for(label i = 0; i < eST.size(); i++)
    {
      
  
      DynamicList<scalar> roadPatchList;

      for(label j=0;j<roads.size();j++)
	{
	  if(emissionGroups[roads[j][1]]==emissionGroups[i])
	     roadPatchList.append(roads[j][0]);	  
	}
      
      
      wordList roadPatchNameList;
      roadPatchNameList.setSize(roadPatchList.size());
      
      
      for(label j=0;j<roadPatchList.size();j++)
	{
	  label roadIndex;
	  roadIndex=label(roads[roadPatchList[j]][0]);
	  if(roadIndex<0 or roadIndex >= roadNames.size())
	     FatalErrorIn(args.executable()) << "Index to roadname given in roadClassification is out of range"<< exit(FatalError);
	  roadPatchNameList[j]=roadNames[roadIndex];
	}

      scalarList roadPatchWeightList;
      roadPatchWeightList.setSize(roadPatchList.size());

      Info<<"debug 5"<<endl;

      for(label j=0;j<roadPatchList.size();j++)
	roadPatchWeightList[j]=roads[roadPatchList[j]][2];      
	
      //word roadListPrefix("patches_");
      //trafficData.lookup(roadListPrefix & emissionGroups[i]) >> roadPatchList;
      
      //word weightListPrefix("weights_");
      //scalarList patchWeightList(trafficData.lookup(weightListPrefix & emissionGroups[i]));
      
      volScalarField& eSTi = eST[i];
      eSTi.internalField()=0;
      scalar totRoadVol;

   
      for(label j = 0; j <roadPatchNameList.size() ; j++)
	{ 
	  Info<<"Processing road patch: "<<roadPatchNameList[j]<<"..."<<endl;
	  //setting eST field for road patch j in emissionGroup i to zero
	  eSTij.internalField()=0.0;
	  
	  //initializing total road volume
	  totRoadVol=0.0;

	  // Declaring patch index labels
	  label patchI;
      
	  // Finding indexes to given patches 
	  patchI=mesh.boundaryMesh().findPatchID(roadPatchNameList[j]);
	  
	  if (patchI == -1)
	    FatalErrorIn(args.executable()) << "Cannot find patch: "<< roadPatchNameList[j] <<" belonging to emissionGroup: "<< emissionGroups[i]<< exit(FatalError);  
  
	  // Declaring and initializing references to boundary patches
	  const polyPatch& pp = mesh.boundaryMesh()[patchI];

	  //Declaring label for cell to be used for information output
	  label infoCellLabel=0;

	  // Looping over all faces of the chosen patch
	  forAll(pp.faceCentres(), patchFaceI)
	    {
	      //get z-coordinate of face centre
	      const point& faceCentre= pp.faceCentres()[patchFaceI];

	      //get label and cell centre z-coordinate to neighbour and owner of the patch face    
	      vector probePoint(faceCentre.x(),faceCentre.y(),faceCentre.z()+0.0001);
            
	      scalar height = emissionHeight + faceCentre.z();
          
	      label nextCell = mesh.findCell(probePoint);
	      scalar cell_z = mesh.C()[nextCell].z();
	      
	      // initialize oldCell variable        
	      label oldCell = nextCell;      
	      
	      //create variable for neighbouring cell centre z-values           
	      vector NbCellCentre;
	      
	      scalar zDist;
	      scalar rDist;
	      scalar angle;

	      label oldNextCell=-9999;
	      infoCellLabel=nextCell;

	      //runs upwards until the cell centre is above the emission height
	      while( cell_z < height)
		{ 
		  if(nextCell==oldNextCell)
		    break;
		  else
		    oldNextCell=nextCell;
		  scalar cellVolume=0.0;
		  cellVolume=mesh.V()[nextCell];
		  // emission source term is set to the roadPatchWeight weight multiplied by
		  //0.001 to get the right fraction of 0.001 kg/s for the patch
		  eSTij.internalField()[nextCell]=roadPatchWeightList[j]*0.001;
		  totRoadVol+=cellVolume;		  
		  //get labels to neighbouring cells 
		  const labelList& cellNeighbours = mesh.cellCells()[nextCell];

		  scalar oldAngle=999;
		  
		  //loop neighbouring cells to find the one most directly above
		  forAll(cellNeighbours, cellI)
		    {
		      NbCellCentre = mesh.C()[cellNeighbours[cellI]];
		      zDist = NbCellCentre.z()-cell_z;
		      rDist = Foam::sqrt(Foam::sqr(NbCellCentre.x())+Foam::sqr(NbCellCentre.y()));
		      angle = Foam::atan(rDist/zDist);
		      
		      if(zDist > 0 && angle < oldAngle )
			{
			  cell_z = NbCellCentre.z();
			  oldAngle = angle;
			  oldCell = nextCell;
			  nextCell = cellNeighbours[cellI];  
			}
		      
		      //dedugging output
		      // Info<< "Face: " << patchFaceI << " cellnb: " << cellI << " cellheight: " <<NbCellCentre.z() << " angle: " << angle << endl;  
		    }
		  
		  if(oldCell == nextCell)      
		    break;
		}
	    }
	  eSTij.internalField()=eSTij.internalField()/totRoadVol;
	  eSTi.internalField()=eSTi.internalField()+eSTij.internalField();
	  Info << "Emissions set for road: " << roadPatchNameList[j]<<endl
	       << tab << "emission group is: " << tab << tab << emissionGroups[i]<<endl
	       << tab << "emission volume is:" << tab << tab <<totRoadVol<<endl
	       << tab << "emission weight is:" << tab << tab <<roadPatchWeightList[j]<<endl
	       << tab << "source term kg/(m3s) is:"<< tab <<eSTi.internalField()[infoCellLabel]<<endl
	       << tab << "corresponding to an emission of : "<<totRoadVol*eSTi.internalField()[infoCellLabel]<<" g/s"<<endl<<endl;
	        
	}
      eST[i].write();	
      }


 Info<< "setEmission finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //
