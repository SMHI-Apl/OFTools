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
  setVehicles

  Description
  Sets vehicle velocity field (Ucar) and vehicle area density (VAD) field

  \*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include "fileName.H"
#include "wallFvPatch.H"
#include "cellSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar lineLength(List<vector> points)
{
  scalar dist=0;
  for(label i=1;i<points.size();i++)
    {
      dist+=Foam::sqrt(mag(points[i]-points[i-1]));
    }
  return dist;
}

//Distance to a line defined by two points
scalar distToLine(vector p1, vector p2, vector p)
{
  scalar t_s;
  vector t_v;
  vector line = p2-p1;
  
  if(p1==p2)
    {
      t_s = 0.0;
    }
  else
    {
      scalar bot = pow(line.x(),2)+pow(line.y(),2)+pow(line.z(),2);
      t_v[0] = (p[0] - p1[0]) * line[0];
      t_v[1] = (p[1] - p1[1]) * line[1];
      t_v[2] = 0;

      t_v[0]=t_v[0]/bot;
      t_v[1]=t_v[1]/bot;
      t_s = max(max(t_v[0],t_v[1]),0.0);
      t_s=min(t_s,1.0);
    }

  line = t_s*line;
  //  line[0]=t_s*line[0];
  //line[1]=t_s*line[1];
  vector pn= p1 + line;
  vector dp = pn-p;
  return mag(dp);
}


//Get direction for line at a given point
vector closestVector(List<vector> points, vector point)
{
  label closestLineIndex=0;
  scalar minDistToLine=9e10;
  for(label i=1;i<points.size();i++)
  {
    scalar dist = distToLine(points[i-1],points[i],point);
    if(dist<minDistToLine)
      {
	closestLineIndex=i-1;
	minDistToLine=dist;
      }
  }    
  vector directionVector = points[closestLineIndex+1] - points[closestLineIndex];
  directionVector = directionVector/mag(directionVector);
  return directionVector;
  
}


int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readTrafficData.H"

  vector UcarVec;
  UcarVec.x()=0;
  UcarVec.y()=0;
  UcarVec.z()=0;

  Ucar.internalField()=UcarVec;

  forAll(roadNames, roadi)
    {
      cellSet roadSet(mesh,roadNames[roadi]+"_vol",IOobject::MUST_READ);
      
      //Calculating the density of the vehicle canopy layer
      //The number of cars trafficing the road at a single moment is set to a constant number
      //The average day traffic (ADT) is used as an estimate together with
      // speed, this gives the momentarily number of cars on the road
      //There is always at least one car on the road when emissions take place, this is set to a minimum number
      //The minima should be weighted with regard to emissions
      //as a simplification, the maximum hourly traffic MHT is a used rather than ADT
      
      scalar speed = roadProperties[roadi][0];
      scalar lanes = roadProperties[roadi][1];
      scalar ADT = roadProperties[roadi][2];
      scalar MHT = 2*ADT/(24.0); //empirical estimate of MHT based on ADT
      scalar roadLength = lineLength(centreLines[roadi]);
      scalar travelTime = roadLength/(speed*1000/3600.0);
      scalar nCarsOnRoad = MHT*travelTime/(3600);      
      nCarsOnRoad = max(1.0, nCarsOnRoad);      
      scalar nUnitLengths = nCarsOnRoad/lanes;
      scalar vol=0;

      
      for(cellSet::iterator cellIter=roadSet.begin();cellIter!=roadSet.end();cellIter++)
	{
	  label i=cellIter.key();
	  vol+=mesh.V()[i];
	}
      
      scalar unitVol = vol/nUnitLengths;
      //dimensionedScalar unitVol("unitVol",(0 3 0 0 0 0 0),unitVolume);
      
      for(cellSet::iterator cellIter=roadSet.begin();cellIter!=roadSet.end();cellIter++)
	{
	  label celli = cellIter.key();
	  vector cellCentre = mesh.C()[celli];
	  vector roadDirection = closestVector(centreLines[roadi],cellCentre);
	  
	  UcarVec.x() = roadDirection.x() * speed*1000/3600.0;
	  UcarVec.y() = roadDirection.y() * speed*1000/3600.0;
	  
	  //dimensionedVector Ucar_init("Ucar_init",(0 1 -1 0 0 0 0),UcarVec);
	  Ucar.internalField()[celli] = UcarVec;
	  VAD.internalField()[celli]=Acar.value()/unitVol;
	  fluidFrac.internalField()[celli]=(unitVol-pow(Acar.value(),3/2))/unitVol;
	}
      Info <<roadNames[roadi]<<"\t\t nrCells "<<roadSet.size()<<"\tspeed "<<speed<<"\tlanes "<<lanes<<"\tMHT "<<MHT<<"\tnCarsOnRoad "<<nCarsOnRoad<<"\troadVol "<<vol<<endl;

    }  
  Ucar.write();
  VAD.write();
  fluidFrac.write();
  Info<< "setVehicle finished successfully\n" << endl;
  return(0);
}
  
  // ************************************************************************* //

