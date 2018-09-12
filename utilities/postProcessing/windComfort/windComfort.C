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
#include <string>
#include <sstream>
#include "PtrList.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "SortableList.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label getDirIndex(SortableList<label> dirs, label dir)
{
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
 #include "postProcess.H"
 #include "setRootCase.H"
 #include "createTime.H"
 #include "createMesh.H"
 #include "createControl.H"
 #include "createFields.H"
 #include "readWindComfortDict.H"
 #include "readWindFields.H"
    
 label nClasses(80); //It is expected that windspeed is always < 40m/s, 0.5m/s intervals for the classes 
 
 scalar classSize(0.5);
 List<scalar> freqTableRow(wdirs.size()+1);//First column is meant for the row sum
 forAll(freqTableRow,coli)
   {
     freqTableRow[coli]=0;
   }

 List< List<scalar> > freqTable(nClasses+1,freqTableRow); //last row is meant for the column sum

 label wspeedClass;
 label wdirClass;
 

 //Count observations for each box and increment frequency table correspondingly
 std::list<tsRow>::iterator row;
 label nCalm;
 label dbg=0;
 for(row=metTS.rows.begin();row!=metTS.rows.end();row++)
   {
     //     Info<<"row"<<dbg++<<endl;
     if(row->data[wdirInd]==double(metTSNodata) or row->data[wspeedInd]<double(detectionLimit))
       nCalm++;
     else

       {
	 wspeedClass= label((scalar(row->data[wspeedInd])-detectionLimit)/classSize);
	 wdirClass = getDirIndex(wdirs,label(row->data[wdirInd]))+1;
	 freqTable[wspeedClass][wdirClass]=freqTable[wspeedClass][wdirClass]+1;
	 freqTable[wspeedClass][0]=freqTable[wspeedClass][0]+1;
	 freqTable[nClasses][wdirClass]=freqTable[nClasses][wdirClass]+1;
       }
   }

 Info<<"Frequency table:"<<endl;
 Info<<"speed\tAll";
 forAll(wdirs,wdiri)
   {
     Info<<"\t"<<wdirs[wdiri];
   }
 Info<<endl;

 forAll(freqTable,rowi)
   {
     if(freqTable[rowi][0]>0)
       Info<<rowi*classSize+detectionLimit<<"-"<<(rowi+1)*classSize+detectionLimit<<"\t";
     forAll(freqTable[rowi],coli)
       {
	 if(freqTable[rowi][0]>0)
	   Info<<freqTable[rowi][coli]<<"\t";
       }
     if(freqTable[rowi][0]>0)
       Info<<endl;
   }	 



 forAll(Udirs,i)
   {
     label refCellIndex;
     refCellIndex = mesh.findCell(refPoints[i]);   
     
     scalar Uref = Udirs[i].internalField()[refCellIndex];
     double perc50;
     scalar minDir,maxDir;
     Info << endl<<"Wind direction: "<<wdirs[i]<<endl;     
     
     if(i==0) {
       minDir=0.5*wdirs[i]-0.5*(360-wdirs[wdirs.size()-1]);
       if(minDir<0)
	 minDir=360+minDir;
     }
     else {
       minDir=wdirs[i]-0.5*(wdirs[i]-wdirs[i-1]);
     }
     if(i==wdirs.size()-1){
       maxDir=wdirs[i]+0.5*(360-wdirs[i]+wdirs[0]);
       if (maxDir>360)
	 maxDir-=360;
     }
     else {
       maxDir=wdirs[i]+0.5*(wdirs[i+1]-wdirs[i]);
     }
 	
     Info<<"Wind direction interval is: "<<minDir<<" < wdir < "<<maxDir<<endl;
     //Calculating median for the current wdir
     perc50 = metTS.condPerc(50, metTS, int(wdirs.size()) , int(wdirInd), int(wspeedInd), double(minDir),double(maxDir),double(metTSNodata));
     scalar fractionWdir(freqTable[nClasses][i+1]/scalar(metTS.nrows));
     Info << "Wind speed 50-percentile:"<<scalar(perc50)<<endl;
     Info << "Percentage of observations for wdir"<<fractionWdir*100<<"%"<<endl;
     Info << "Uref: "<<refPoints[i]<<" is: "<<Uref<<endl;
     if(Uref<0.5)
       Info<<"!Warning, Uref is close to zero, check coordinates of reference-points"<<endl;

       
     volScalarField Uscaled = (scalar(perc50)/Uref) * Udirs[i] + Umin;
     Info<<"Scaled U in ref point: "<<Uscaled.internalField()[refCellIndex]<<endl;
     volScalarField kscaled=scalar(perc50)/Uref * kdirs[i];
     Info<<"Scaled k in ref point: "<<kscaled.internalField()[refCellIndex]<<endl;
     volScalarField I = sqrt(2/3.0*kscaled)/Uscaled;
     Info<<"Turbulent intensity, I, in ref point: "<<I.internalField()[refCellIndex]<<endl;
     volScalarField Uekv = Uscaled*(1+alpha*I)/(1+alpha*Iref);
     Info<<"Uekv for median of current wdir in ref point: "<<Uekv.internalField()[refCellIndex]<<endl;
     UekvMedian = UekvMedian + fractionWdir*Uekv; 
     Info<<"Accumulated UekvMedian in ref point: "<<UekvMedian.internalField()[refCellIndex]<<endl;
     

     for(label j=0;j<nClasses;j++)
       {
	 scalar wspd=j*classSize+classSize*0.5+detectionLimit;
	 if(freqTable[j][0]>0) //If there are any observations within the current wspeed class
	   {	
	     //Info<<"Wind speed: "<<wspd<<" m/s"<<", found "<<freqTable[j][0]<<" observations"<<endl;
	     Uscaled=wspd/Uref * Udirs[i]+Umin;
             //Info<<"Scaled U for ref wspd "<<wspd<<" in ref point: "<<Uscaled.internalField()[refCellIndex]<<endl;
             kscaled=wspd/Uref * kdirs[i];
             //Info<<"Scaled k for ref wspd "<<wspd<<" in ref point: "<<kscaled.internalField()[refCellIndex]<<endl;
	     I=sqrt(2/3.0*kscaled)/Uscaled;
             //Info<<"Turbulent intensity, I, for ref wspd "<<wspd<<" in ref point: "<<I.internalField()[refCellIndex]<<endl;
	     Uekv=Uscaled*(1+alpha*I)/(1+alpha*Iref);
             //Info<<"Uekv for ref wspd "<<wspd<<" in ref point: "<<Uekv.internalField()[refCellIndex]<<endl;
	     
	     forAll(Uekv.internalField(),celli)
	       {
		 if(Uekv.internalField()[celli] >= criteriaUlim)
		   UekvProc.primitiveFieldRef()[celli] = UekvProc.internalField()[celli] + freqTable[j][i+1];
	       }      
	   }
       }
     Info<<"UekvProc in ref point: "<<UekvProc.internalField()[refCellIndex]/scalar(metTS.nrows)*100<<"%"<<endl;
   }

 UekvProc = UekvProc / scalar(metTS.nrows)*100;
 UekvProc.write();
 UekvMedian.write();
 Info<< "windComfort finished successfully\n" << endl;
 return(0);
}

// ************************************************************************* //
