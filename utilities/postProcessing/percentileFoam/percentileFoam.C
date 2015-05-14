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
#include "ConcFileArchive.H"
#include "tsRow.H"
#include "timeSeries.H"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "PtrList.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readStatisticalData.H"
#   include "readEmissions.H"  
#   include "createFields.H"
#   include "createTimeSeries.H"
  // * ** * * * * * * * * * * * * * * * * * * * * * * //  
  
  //Loops over all hours, calculating a weight corresponding to the current meteorology and emissions
  //The weights are stored in a timeseries with one column per file in the ConcFileArchive 
  Info<<"Loops over all hours to calulate weights for dispersion cases"<<endl;
  emisHour=emisTS.rows.begin();
  for(hour=metTS.rows.begin();hour!=metTS.rows.end();hour++)
    {
      //Gets the vector with weights caused by the meteorology
      //(linear interpolation between windspeeds, gaussian wdir variation with std sigma during single hours)
      std::vector<double> metWeights = ca.getConcMetWeights((*hour).data[wspeedInd],(*hour).data[wdirInd],(*hour).data[sigmaInd]);
      std::vector<double> emisWeights(metWeights.size(),0);
      std::vector<double> concNameWeights;
      
      //Iterates over all concNames in the ConcFileArchive
      //gets a weight vector with a weight for each ConcFileArchive (weight=1 if the file refers to the current concName, otherwise 0)
      //For each concName the emission is also extracted from the emission timeseries and multiplied with the weight vector to give an emission weight vector
      //The resulting emission weights are summed over all concNames to get a complete emission weight vector
      //Finally the emission weight vector is multiplied with the meteorological weight vector to give a complete weight vector
      //The emission timeseries is presumed to be in units mg/s (a conversion is made to match the unit emission of 1 g/(s sourceGroup) on which the concentration fields are based)

      std::vector<std::string>::iterator concName;
      for(concName=ca.concNames.begin();concName!=ca.concNames.end();concName++)
        {
          concNameWeights=ca.getConcNameWeights(*concName);
          int concNameInd=emisTS.getColInd(*concName);
          if(concNameInd==-1)
            {
              Info<<"Error: parameter "<<*concName<<" not found in timeseries header"<<endl;
              FatalErrorIn(args.executable())
                << "Error: parameter "<< *concName<<" not found in timeseries header"<< exit(FatalError);
            }
          emisWeights=add(emisWeights,multi(concNameWeights,(*emisHour).data[concNameInd]/1000.0));
        }
      tsRow hourWeights((*hour).year, (*hour).month, (*hour).day, (*hour).hour,multi(metWeights,emisWeights));  
      weightTS.appendRow(hourWeights);

      tsRow debugWeightRow((*hour).year, (*hour).month, (*hour).day, (*hour).hour,metWeights);
      debugTS.appendRow(debugWeightRow);
      emisHour++;
      Info<<"emissions are: "<<emisHour->year<<"-"<<emisHour->month<<"-"<<emisHour->day<<" "<<emisHour->hour<<endl;

  }

  weightFid<<debugTS;
  weightFid.close();

  volScalarField& hourlyAvg=statFields[0];
  
  Info<<"Preprocessing of met-data and emissions ready!"<<endl;
  Info<<"Starting loop over cells!..."<<endl;

  //Factor to convert from massfraction to ug/m3
  //temporary fix with hard-coded air density
  double typeConversionFactor=1.225*1.0e9;


  //Creating quality control files
  std::ofstream qcFid;
  fileName qcFileName=args.rootPath()+"/"+args.globalCaseName()+"/"+"qualityControlFile.txt";
  
  Info<<"quality control file: "<<qcFileName.c_str()<<endl;
  
  qcFid.open(qcFileName.c_str());
  if(!qcFid.is_open())
    FatalErrorIn(args.executable())
      << "Could not open quality control file for max conc"<< exit(FatalError);
  qcFid<<"Quality assurance of calculations\n";

  //Initializing progress counters
  scalar cellCounter=1;
  label nCells=hourlyAvg.internalField().size();
  scalar doneProc1=0;
  scalar doneProc2=0;
  int probesDone=0;
  
  //Iterating over all internal cells
  Info<<" Number of cells is: "<<nCells<<endl;
  forAll(hourlyAvg.internalField(),celli)
    {
      //Progress counter
      doneProc2=cellCounter/scalar(nCells)*100.0;
      if(doneProc2-doneProc1>1)	
	{
	  Info<<" Done "<<int(doneProc2)<<" %"<<endl;
	  doneProc1=doneProc2;
	}

      //Getting cell values from archived conc fields
      ca.nextCellValue();

      //Checking if cell label match any of the probe cell labels
      int foundProbe=0;
      int probeNr;
      for(probeNr=0;probeNr<int(probeCells.size());probeNr++)
	{
	  if(probeCells[probeNr]==celli)
	    {
	      Info<<"Found probe cell, writing timeseries to extractedTimeseries.txt"<<endl;
	      foundProbe=1;
	      probesDone++;
	      break;
	    }
	}

      //Initializing average calculations
      double sumConc=0.0;
      double nHours=0.0;
   
      //Start timeseries calculations for current cell
      if(write3DFields or foundProbe)
	{
	  bgHour=bgTS.rows.begin();
	  concHour=concTS.rows.begin();
	  if(foundProbe)
	    probeHour=probeTS.rows.begin();
	  for(hour=weightTS.rows.begin();hour!=weightTS.rows.end();hour++)
	    {
	      (*concHour).data[0]=sum(multi((*hour).data,ca.concs))*typeConversionFactor;	      
	      //Checking that calculated concentration is not greater than maxConc
	      //if so writing individual concentrations to file for closer checking
	      if( (*concHour).data[0]>maxConc or ( foundProbe and hour==weightTS.rows.begin() ) )
		{
		  std::vector<ScalarFoamFile>::iterator concFile;
		 
		  if(foundProbe)
		    {
		      Info<<"Quality control data written to file: "<<qcFileName<< endl;
		      qcFid<<"Probe name: "<<probeNames[probeNr]<<"\n";
		    }
		  else
		    Info<<"Warning: "<<"There are unvalid concentration fields, sum of conc in cell "<<cellCounter<<" is: "<< (*concHour).data[0] <<", max allowed conc is set to: "<<maxConc<<"summary written to: "<<qcFileName<< endl;
		  
		  qcFid<<"Time is: "<<hour->year<<"-"<<hour->month<<"-"<<hour->day<<"-"<<hour->hour<<"\n";
		  qcFid<<"Filename\torgConc\tweight\ttypeConversion\ttotConc\n";
		  int concInd=0;
		  for(concFile=ca.fileVec.begin();concFile!=ca.fileVec.end();concFile++)
		    {
		      qcFid<<concFile->getDirName()<<" "<<concFile->fileName<<"\t";
		      qcFid<<ca.concs[concInd]<<"\t";
		      qcFid<<(*hour).data[concInd]<<"\t";
		      qcFid<<typeConversionFactor<<"\t";
		      qcFid<<ca.concs[concInd]*(*hour).data[concInd]*typeConversionFactor<<"\n";		  
		      concInd++;
		    }
		}
	      

	      (*concHour).year=(*hour).year;
	      (*concHour).month=(*hour).month;
	      (*concHour).day=(*hour).day;
	      (*concHour).hour=(*hour).hour;

	      if(background)
		(*concHour).data[0]+=(*bgHour).data[0];
	      
	      sumConc+=(*concHour).data[0];
	      nHours+=1.0;

	      if(foundProbe)
		{
		  (*probeHour).data[probeNr]=(*concHour).data[0];
		  probeHour++;
		}

	      bgHour++;
	      concHour++;
	    }
	  cellCounter++;
	  //std::ofstream testOut2;
	  //testOut2.open("/local_disk/dsegerss/OpenFOAM/dsegerss-1.3/run/slbProject/Hornsgatan_test/testOut2.txt");
	
	  if(probesDone==int(probeCells.size()) and not write3DFields)
	    break;          
	}
      
      if(write3DFields)
	{
	  if(indDailyAvgPerc)
	    {
	      timeSeries dailyAvgTS;
	      dailyAvgTS=concTS.dailyAverage(int(nValidsForDaily));
	      dailyAvgTS.header.push_back("dailyAverage");
	      dailyAvgTS.sort();      
	      statFields[indDailyAvgPerc].internalField()[celli]=dailyAvgTS.percentile(dailyAvgPercentile);
	      
	    }

	  if(indDailyMaxPerc)
	    {
	      
	      timeSeries dailyMaxTS;
	      dailyMaxTS=concTS.dailyMax(int(nValidsForDaily));
	      dailyMaxTS.header.push_back("dailyMax");
	      dailyMaxTS.sort();
	      statFields[indDailyMaxPerc].internalField()[celli]=dailyMaxTS.percentile(dailyMaxPercentile);
	    }

	  if(indHourlyPerc)
	    {
	      concTS.sort();
	      //statFields[indHourlyPerc].internalField()[celli]=concTS.percentile(hourlyPercentile);;
	     
	      //Debugging
	      scalar cellHourlyPerc=concTS.percentile(hourlyPercentile);
	      statFields[indHourlyPerc].internalField()[celli]=cellHourlyPerc;
	      //if(cellHourlyPerc>5.71e5)
	      //	{
	      //	  Info<<"Cell index for hourly percentile > 5.71e5 is: "<< cellCounter-1<<endl;
	      //  int last_ind=concTS.nrows-1;
	      //  Info<<"Date is: "<<concTS.getYear(last_ind)<<"-"<<concTS.getMonth(last_ind)<<"-"<<concTS.getDay(last_ind)<<" "<<concTS.getHour(last_ind)<<" value:"<< concTS.getValue(last_ind)<<endl;
	      //	}
	    }
  

	  hourlyAvg.internalField()[celli]=sumConc/nHours;
	}
      cellCounter++;
    }

  qcFid.close();

  if(write3DFields)
    {
      //Writing results
      //      hourlyAvg=hourlyAvg;
      hourlyAvg.write();
      
      if(indHourlyPerc)
	statFields[indHourlyPerc].write();
      
      if(indDailyAvgPerc)
	statFields[indDailyAvgPerc].write();


      if(indDailyMaxPerc)
	statFields[indDailyMaxPerc].write();
    }

      if(int(probeCells.size()>0))
	{	  
	  probeFid<< probeTS;
	  probeFid.close();
	}

  Info<< "percentileFoam finished successfully\n" << endl;
  return(0);
}

// ************************************************************************* //
