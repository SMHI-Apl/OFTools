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
  percentileFoam

  Description
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include "tsRow.H"
#include "timeSeries.H"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "PtrList.H"
#include <iostream>
#include "Raster.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "SortableList.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label getDirIndex(SortableList<label> dirs, label dir)
{
  dirs.sort();
  label minDist=360;
  label dirIndex=-1;
  forAll(dirs,diri)
    {
      scalar dist=min(abs(dirs[diri]-dir),abs(abs(dirs[diri]-dir)-360));
      if( dist<minDist)
	{
	  minDist=dist;
	  dirIndex=diri; 
	}
    }
  return dirIndex;
}
  

int main(int argc, char *argv[])
{  
#   include "setRootCase.H"
#   include "createTime.H"
#   include "readWindStatisticsDict.H"

  PtrList<Raster> results(6);
  forAll(results,resi)
    {
      results.set(resi,new Raster(rasters[0],0));
      results[resi]=results[resi]+(rasters[0]*0);
    }
  
  forAll(rasters,wdiri)
    {
      Info<<"Raster max before scaling: "<<rasters[wdiri].max();
      rasters[wdiri].scale(1/refWinds[wdiri]);
      Info<<", after scaling: "<<rasters[wdiri].max()<<endl;
    }
  Info<<"Normalized wind rasters with reference winds"<<endl;



  //Extracting time series in probe points
  std::ofstream probeFid;
  if(probePointList.size()>0)
    {

      Info<<"probe file: "<<probeFileName.c_str()<<endl;
      probeFid.open(probeFileName.c_str());
      if(!probeFid.is_open())
        FatalErrorIn(args.executable())
          << "Could not open probe output file"<< exit(FatalError);
    }

  timeSeries probeTS;  
  //Creating an empty ts row to initialize the probe timeseries with
  tsRow probeRow(int(probeNameList.size())+4);
  for(int i=0;i<probeRow.ncols;i++)
    probeRow.data[i]=-9999.0;
  probeTS.header.push_back("wdir_coast");
  probeTS.header.push_back("wspeed_coast");
  probeTS.header.push_back("wspeed_vanern");
  probeTS.header.push_back("wdir_vanern");
  std::list<tsRow>::iterator hourMesanLine,hourWaspVanern,hourWaspCoast,hourMesanVanern,hourProbe;

  //Filling the probe time series and the concentration time series with dates from the meteorology time series
  for(hourMesanLine=mesanLineTS.rows.begin();hourMesanLine!=mesanLineTS.rows.end();hourMesanLine++)
    { 
      probeRow.year=(*hourMesanLine).year;
      probeRow.month=(*hourMesanLine).month;
      probeRow.day=(*hourMesanLine).day;
      probeRow.hour=(*hourMesanLine).hour;
      probeTS.appendRow(probeRow);
    }

  forAll(probePointList,probei)
    {
      probeTS.header.push_back(std::string(probeNameList[probei].c_str()));
      List<scalar> extrU(rasters.size());
      scalar extrUval;
      forAll(rasters,rasti)
	{
	  extrU[rasti]=rasters[rasti].getValue(probePointList[probei].x()-subtractedX,probePointList[probei].y()-subtractedY);
	  Info<<"U for res nr "<<rasti<<" is: "<<extrU[rasti]<<endl;
	}
      
      hourMesanVanern=mesanVanernTS.rows.begin();
      hourWaspCoast=waspCoastTS.rows.begin();
      hourProbe=probeTS.rows.begin();

      for(hourMesanLine=mesanLineTS.rows.begin();hourMesanLine!=mesanLineTS.rows.end();hourMesanLine++)
	{
	  label wdir = label(hourMesanVanern->data[mesanVanernWdirInd]);
	  scalar Ulocal = scalar(hourMesanVanern->data[mesanVanernWspeedInd]);
	  scalar Ucoast =scalar(hourWaspCoast->data[waspCoastWspeed100Ind]);
	  scalar wdirCoast=scalar(hourWaspCoast->data[waspCoastWdirInd]);
	  label resInd=-1;
	  label rasterIndex = getDirIndex(wdirs,wdir);
	  extrUval = extrU[rasterIndex]*Ulocal;	  		  
	  hourProbe->data[0]=wdirCoast;
	  hourProbe->data[1]=Ucoast;
	  hourProbe->data[2]=Ulocal;
	  hourProbe->data[3]=wdir;
	  hourProbe->data[probei+4]=extrUval;
	
	  hourWaspCoast++;
	  hourMesanVanern++;
	  hourProbe++;
	}
    }

  probeFid<<probeTS;
  probeFid.close();
 
  hourMesanVanern=mesanVanernTS.rows.begin();
  hourWaspCoast=waspCoastTS.rows.begin();
  hourWaspVanern=waspVanernTS.rows.begin();
 
  for(hourMesanLine=mesanLineTS.rows.begin();hourMesanLine!=mesanLineTS.rows.end();hourMesanLine++)
   {
     label wdir = label(hourMesanVanern->data[mesanVanernWdirInd]);
      scalar Ulocal = scalar(hourMesanVanern->data[mesanVanernWspeedInd]);
      scalar Ucoast =scalar(hourWaspCoast->data[waspCoastWspeed100Ind]);
      scalar clouds = hourMesanLine->data[cloudInd];
      scalar T = hourMesanLine->data[TInd];
      label resInd=-1;

      if(T>=20)
	{
	  if( clouds <=15 and Ucoast>=13)
	    {
	      Info<<"time: "<<hourMesanLine->year<<"-"<<hourMesanLine->month<<"-"<<hourMesanLine->day<<"-"<<hourMesanLine->hour<<"\t"
		  <<"wdir: "<<wdir<<"\t"
		  <<"Ucoast: "<<Ucoast<<"\t"
		  <<"Ulocal: "<<Ulocal<<"\t"
		  <<"T: "<<T<<"\t"
		  <<"cloud: "<<clouds<<"\t";
	      label rasterIndex = getDirIndex(wdirs,wdir);
	      Raster hourRes;
	      hourRes=rasters[rasterIndex]*double(Ulocal);
	      Info<<"rastMax: "<<hourRes.max()<<"!"<<endl; 
	      
	      if(T>=20)
		{
		  results[0].whereAdd( hourRes<=0.6 ,1 );
		  results[1].whereAdd( hourRes<=1.5 ,1 );
		  results[2].whereAdd( hourRes<=3.0 ,1 );
		}
	      
	      if(T>= 25)
		{
		  results[3].whereAdd( hourRes<=0.6 ,1 );
		  results[4].whereAdd( hourRes<=1.5 ,1 );
		  results[5].whereAdd( hourRes<=3.0 ,1 );
		}
	    }
	}
    
      hourWaspVanern++;
      hourWaspCoast++;
      hourMesanVanern++;
    }
   

  forAll(results,resi)
    {
      OStringStream resFileName, sampleFileName;
      resFileName<<"case_"<<resi<<".asc";
      sampleFileName<<"case_"<<resi<<"_samples.asc";
      fileName resFile=resultDir/resFileName.str();
      fileName sampleFile=resultDir/sampleFileName.str(); 
      OFstream pointSampleFile(sampleFile);
      
      results[resi].xll+=subtractedX;
      results[resi].yll+=subtractedY;
      results[resi].xur+=subtractedX;
      results[resi].yur+=subtractedY;

      results[resi].write(resFile.c_str());
      Info<<"Wrote result "<<resi<<" to disk"<<endl;
      Info<<"Extracting points along line"<<endl;
      pointSampleFile << "X\tY\tdistance\tobservations\n";

      scalar stepsize=50;
      scalar dist=0;
      scalar accDist=0;
      for(label i=1;i<points.size();i++)
	{
	  Info<<"sampling linepart: "<<i<<endl;
	  vector p2=points[i];  
	  vector p1=points[i-1];
	  vector lineVec=p2-p1;
	  scalar length=mag(lineVec);
	  vector dirVec=lineVec/length;
	  label steps = label(length/scalar(stepsize)+1);
	  for(label j=0;j<steps;j++)
	    {
	      vector samplePoint=p1+dirVec*(j*stepsize);
	      dist=accDist+(mag(samplePoint-p1));
	      pointSampleFile<<samplePoint.x()<<"\t"<<samplePoint.y()
			     <<"\t"<<dist<<"\t"
			     <<results[resi].getValue(samplePoint.x(),samplePoint.y())<<"\n";
	      // Info<<"samplePoint: "<<samplePoint.x()<<"\t"<<samplePoint.y()
	      //	     <<"\t"<<dist<<"\t"
	      //		     <<results[resi].getValue(samplePoint.x(),samplePoint.y())<<"\n";
	    }
	  accDist+=length;
	}
      //      pointSampleFile.close();
    }

  

  Info<< "windStatistics finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //
