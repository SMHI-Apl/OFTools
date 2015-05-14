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
#include "cellSet.H"

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

  //Looping over all species(emissionGroups)
  for(label i = 0; i < emissionGroups.size(); i++)
    {

      Info << "Processing emission group: " << emissionGroupNames[i] << endl;

      DynamicList<label> roadSetList;
      DynamicList<word> roadSetNameList;
      DynamicList<scalar> roadSetWeightList;

      for(label j=0;j<roads.size();j++)
	{
	  label roadEmissionGroup=label(roads[j][1]);
	  label roadIndex=label(roads[j][0]);
	  
	  if(roadIndex<0 or roadIndex >= roadNames.size())
	    FatalErrorIn(args.executable()) << "Index to roadname given in roadClassification is out of range"<< exit(FatalError);
	  
	  scalar roadWeight=roads[j][2];
	  word roadName=roadNames[roadIndex];
	  
	  
	  if(roadEmissionGroup==i and roadWeight>0.0 )
	    {
	      roadSetList.append(roadIndex);
	      roadSetNameList.append(roadName);
	      roadSetWeightList.append(roadWeight);
	    }
	}

      Info << endl <<"Emission group summary for: " << emissionGroupNames[i] << endl
	   << "Roads with non zero weights: "<<endl;
      scalar weightSum=0;
      for(label j=0;j<roadSetNameList.size();j++)
	{
	  Info << tab << j+2 << ": " << roadSetNameList[j]<< " weight: " << roadSetWeightList[j] << endl;
	  weightSum=weightSum+roadSetWeightList[j];
	}
      Info << "Total weight is: " << weightSum << endl << endl;
      
      //eSTi is a reference to the current emissionGroup source term
      volScalarField& eSTi = eST[i];
      eSTi.internalField()=0;
      scalar totRoadVol;
   
      for(label j = 0; j <roadSetNameList.size() ; j++)
	{ 
	  Info<<"Processing road set: "<<roadSetNameList[j]<<"..."<<endl;
	  //setting eST field for road patch j in emissionGroup i to zero
	  eSTij.internalField()=0.0;
	  
	  //initializing total road volume
	  totRoadVol=0.0;

	  cellSet roadSet(mesh,roadSetNameList[j]+word("_vol"),IOobject::MUST_READ);
	  label celli;

	  //Looping over all cells in roadSet
	  for(cellSet::iterator cellIter=roadSet.begin();cellIter!=roadSet.end();cellIter++)
	    {
	      celli = cellIter.key();
	      // emission source term is set to the roadSetWeight weight multiplied by
	      //0.001 to get the right fraction of 0.001 kg/s for the set
	      eSTij.internalField()[celli]=roadSetWeightList[j]*0.001;
	      totRoadVol+=mesh.V()[celli];
	    }

	  //Dividing with volume to get emission in units kg/(m3*s)
	  eSTij=eSTij/totRoadVol;

	  //adding contribution for road to the emissionGroup
	  eSTi=eSTi+eSTij;

	  Info << "Emissions set for road: " << roadSetNameList[j]<<endl
	       << tab << "emission group is: " << tab << tab << emissionGroups[i]<<endl
	       << tab << "emission volume is:" << tab << tab <<totRoadVol<<endl
	       << tab << "emission weight is:" << tab << tab <<roadSetWeightList[j]<<endl
	       << tab << "source term kg/(m3s) in cell nr "<<celli<<" is:"<< tab <<eSTij.internalField()[celli]<<endl
	       << tab << "corresponding to an emission of : "<<1e6*totRoadVol*eSTij.internalField()[celli]<<" mg/s"<<endl<<endl;
	}
      scalar maxEmission=0;
      label nvalues=0;
      
      forAll(eSTi.internalField() ,cell_i)
	{
	  scalar val=eSTi.internalField()[cell_i];
	  if(val>0)
	    nvalues++;
	      
	  if(val>maxEmission)
	    maxEmission=val;
	}   
		
      Info<< "Max emission for "<<emissionGroupNames[i]<<": "<<maxEmission<< endl;
      Info<< "Number of cells with value > 0 :"<<nvalues<<endl;

      eSTi.write();	
      }


 Info<< "setEmission finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //
