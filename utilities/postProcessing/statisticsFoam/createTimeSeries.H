  ConcFileArchive ca(concFileNames);  
  timeSeries metTS(metFile.c_str());
  timeSeries emisTS(emisFile.c_str());    
  timeSeries bgTS;
  if(background)
      bgTS = timeSeries(bgFile);
  
  Info<<"Emission hours read from input file: "<<emisTS.nrows<<endl;
  emisTS.removeNodata();
  Info<<"Emission hours with complete data: "<<emisTS.nrows<<endl;
  Info<<"Meteorology hours read from input file: "<<metTS.nrows<<endl;
  metTS.removeNodata();
  Info<<"Meteorology hours with complete data: "<<metTS.nrows<<endl;

  if(background)
    {
      Info<<"Background hours read from input file: "<<bgTS.nrows<<endl;
      bgTS.removeNodata();
      Info<<"Background hours with complete data"<<bgTS.nrows<<endl;

      //Match the dates of all time series (make sure all dates exist for all timeSeries)
      metTS.matchDates(bgTS);
      bgTS.matchDates(metTS);
      bgTS.matchDates(emisTS);
    }
  emisTS.matchDates(metTS);
  metTS.matchDates(emisTS);
  
  Info<<"Hours to be used for statistics: " << emisTS.nrows<<endl;
  
  if(metTS.nrows != emisTS.nrows)
    {
      FatalErrorIn(args.executable())
        << "The number of rows in the meteorological timeseries ("
	<< metTS.nrows
	<< ") differs from the number of rows in the emission timeseries ("
	<< emisTS.nrows<<")"
	<< exit(FatalError);
    }

    if(background && bgTS.nrows != emisTS.nrows)
    {
      FatalErrorIn(args.executable())
        << "The number of rows in the background timeseries differs "
	<< "from the number of rows in the emission timeseries and "
	<< "the meteorological timeseries"
	<< exit(FatalError);
    }

  timeSeries weightTS;
  std::vector<ScalarFoamFile>::iterator fid;
  for(fid=ca.fileVec.begin();fid!=ca.fileVec.end();fid++)
    weightTS.header.push_back((*fid).fileName);
  
  //Finding probe cells and creating time series for probe points  
  //Creating an empty probe time series row to initialize the probe time series with 
  tsRow probeRow(int(probeList.size()));
  for(int i=0;i<probeRow.ncols;i++)
    probeRow.data[i]=-9999.0;

  //creating an empty concentration time series row to initialize the conentration timeseries
  tsRow concRow(1);
  concRow.data[0]=-9999.0;

  //Creating a vector of probe cell labels and a probe time series, setting the time series headers to the probe names
  timeSeries probeTS;
  std::vector<label> probeCells;  
  forAll(probeList,probe)
    {
      label probeCellLabel;
      probeCellLabel = mesh.findCell(probeList[probe]);
      if(probeCellLabel<0)
	FatalErrorIn(args.executable())
	  << "The probe named: "<<probeNames[probe]<<" is outside of mesh" << exit(FatalError);
      probeCells.push_back(probeCellLabel);
      Info<<"Found probe cell: "<< probeNames[probe]<<" at index: "<<probeCellLabel<<endl;
      probeTS.header.push_back(std::string(probeNames[probe].c_str()));
    }

  std::ofstream probeFid;
  if(int(probeCells.size()>0))
    {
      Info<<"probe file: "<<probeFile.c_str()<<endl;
      probeFid.open(probeFile.c_str());
      if(!probeFid.is_open())
        FatalErrorIn(args.executable())
          << "Could not open probe output file"<< exit(FatalError);
    }

  //Creating a initialized concentration time series
  timeSeries concTS;
  concTS.header.push_back("concentration[g/m3]");

  //Filling the probe time series and the concentration time series with dates from the meteorology time series
  std::list<tsRow>::iterator rowIter;
  for(rowIter=metTS.rows.begin();rowIter!=metTS.rows.end();rowIter++)
    { 
      probeRow.year=(*rowIter).year;
      probeRow.month=(*rowIter).month;
      probeRow.day=(*rowIter).day;
      probeRow.hour=(*rowIter).hour;
      probeTS.appendRow(probeRow);

      concRow.year=(*rowIter).year;
      concRow.month=(*rowIter).month;
      concRow.day=(*rowIter).day;
      concRow.hour=(*rowIter).hour;
      concTS.appendRow(concRow);
    }

  //preprocesing weights to get a time series with weights
  //for all concentration files all the hours
  std::list<tsRow>::iterator hour, emisHour,bgHour,concHour,probeHour, metHour;
  int wdirInd=metTS.getColInd("wdir");
  int wspeedInd=metTS.getColInd("wspeed");
  int sigmaInd=metTS.getColInd("sigmawdir");
  int globInd = metTS.getColInd("glob");
  int tempInd = metTS.getColInd("temp");

  //Gunnars metod: sigmaInd=min(0.5/u*180/pi,19.3*pow(u,-0.327));
  if(wdirInd==-1)
    {
      Info<<"Error: parameter wdir not found in timeseries header"<<endl;
      FatalErrorIn(args.executable())
        << "Error: parameter wdir not found in timeseries header"<< exit(FatalError);
    }
  if(wspeedInd==-1)
    {
      Info<<"Error: parameter wspeed not found in timeseries header"<<endl;
      FatalErrorIn(args.executable())
        << "Error: parameter wspeed not found in timeseries header"<< exit(FatalError);
    }
  if(sigmaInd==-1)
    {
      Info<<"Error: parameter sigmawdir not found in timeseries header"<<endl;
      FatalErrorIn(args.executable())
        << "Error: parameter sigmawdir not found in timeseries header"<< exit(FatalError);
    }
  
if (NO2_conversion)
  {
    if(globInd==-1)
      {
	Info<<"Error: parameter glob (required for NOx to NO2 conversion),"
	    <<" not found in timeseries header"<<endl;
	FatalErrorIn(args.executable())
	  << "Error: parameter glob (required for NOx to NO2 conversion)"
	  <<"not found in timeseries header"<< exit(FatalError);
      }
    if(tempInd==-1)
      {
	Info<<"Error: parameter temp (required for NOx to NO2 conversion),"
	    <<" not found in timeseries header"<<endl;
	FatalErrorIn(args.executable())
	  << "Error: parameter temp (required for NOx to NO2 conversion)"
	  <<"not found in timeseries header"<< exit(FatalError);
      }
  }

int NOxbInd = bgTS.getColInd("NOx_b");
int NO2bInd = bgTS.getColInd("NO2_b");
int O3bInd = bgTS.getColInd("O3_b");

if(NO2_conversion)
  {
    if(NOxbInd == -1)
      {
	Info<<"Error: parameter NOx_b (required for NOx to NO2 conversion),"
	    <<" not found in background timeseries header"<<endl;
	FatalErrorIn(args.executable())
	  << "Error: parameter NOx_b (required for NOx to NO2 conversion)"
	  <<"not found in background timeseries header"<< exit(FatalError);
      }

    if(NO2bInd == -1)
      {
	Info<<"Error: parameter NO2_b (required for NOx to NO2 conversion),"
	    <<" not found in background timeseries header"<<endl;
	FatalErrorIn(args.executable())
	  << "Error: parameter NO2_b (required for NOx to NO2 conversion)"
	  <<"not found in background timeseries header"<< exit(FatalError);
      }

    if(O3bInd == -1)
      {
	Info<<"Error: parameter O3_b (required for NOx to NO2 conversion),"
	    <<" not found in background timeseries header"<<endl;
	FatalErrorIn(args.executable())
	  << "Error: parameter O3_b (required for NOx to NO2 conversion)"
	  <<"not found in background timeseries header"<< exit(FatalError);
      }
  }
   



  //Creating debug output file for weights
  std::ofstream weightFid;
fileName weightFile=args.rootPath()+"/"+args.globalCaseName()+"/"+"weightDebugFile.txt";

  Info<<"weight debug file: "<<weightFile.c_str()<<endl;
  weightFid.open(weightFile.c_str());
  if(!weightFid.is_open())
    FatalErrorIn(args.executable())
      << "Could not open weight debug file"<< exit(FatalError);

  timeSeries debugTS(weightTS);

