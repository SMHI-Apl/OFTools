  for(int i=0;i<int(concFileNames.size()); i++)
      Info<<concFileNames[i]<<endl;

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

  Info<<"Included "<<int(ca.fileVec.size())<<" files from concArchive"<<endl;
  if(ca.fileVec.size() == 0) {
      FatalErrorIn(args.executable())
          << "No conc files included "
          << "(check format of file list in statisticsDict)"         
          << exit(FatalError);
  }
  //Finding probe cells and creating time series for probe points  
  //Creating an empty probe time series row to initialize the probe time series with 
  tsRow probeRow(1);  
  probeRow.data[0]=-9999.0;

  //creating an empty concentration time series row to initialize the conentration timeseries
  tsRow concRow(1);
  concRow.data[0]=-9999.0;

  //Creating a vector of probe cell labels and a probe time series, setting the time series headers to the probe names
  Info<<"Locating probe points in mesh"<<endl;
  PtrList<timeSeries> probeTsList;
  std::vector<label> probeCells;  
  forAll(probeList,probe)
  {
      label probeCellLabel;
      probeCellLabel = mesh.findCell(probeList[probe]);
      if(probeCellLabel >=0 ) {
          probeCells.push_back(probeCellLabel);
          std::cout <<"Found probe cell: "<< probeNames[probe]<<" at index: "<<probeCellLabel<<"\n";
          probeTsList.append(new timeSeries());
          (probeTsList.last()).header.push_back(std::string(probeNames[probe].c_str()));
      }
      //FatalErrorIn(args.executable())
      //<< "The probe named: "<<probeNames[probe]<<" is outside of mesh" << exit(FatalError);
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
      for(PtrList<timeSeries>::iterator probeIter=probeTsList.begin();
          probeIter!=probeTsList.end();probeIter++) {
          (*probeIter).appendRow(probeRow);
      }

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
   

  //Creating debug output file for weight
  std::ofstream weightFid;
  fileName weightFile=args.rootPath()+"/"+args.globalCaseName()+"/"+"weightDebugFile.txt";
  Info<<"weight debug file: "<<weightFile.c_str()<<endl;
  timeSeries debugTS(weightTS);

  //Creating quality control files
  std::ofstream qcFid;
  fileName qcFileName=args.rootPath()+"/"+args.globalCaseName()+"/"+"qualityControlFile.txt";  
  Info<<"quality control file: "<<qcFileName.c_str()<<endl;

  std::ofstream ovFid;
  fileName overviewFileName=args.rootPath()+"/"+args.globalCaseName()+"/"+"concFieldOverview.csv";
  Info<<"Conc field overview file: "<<overviewFileName.c_str()<<endl;
  // Clear old output file that will be appeded by parallel writes

  // Clear old output files that will be appeded by parallel writes
  if( Pstream::master()) {
    ovFid.open(overviewFileName.c_str());
    if(!ovFid.is_open())
      FatalErrorIn(args.executable())
	<< "Could not open file " << overviewFileName << "for writing"
	<< exit(FatalError);
    
    ovFid<<"# Concentration field overview for case: "<< args.globalCaseName()<<"\n";
    // write field names to header in order of archive fileVec
    ovFid<<"Probe-name";
    for(int concFieldInd=0;concFieldInd < int(ca.fileNames.size()); concFieldInd++) {
        fileName fieldFilePath(ca.fileNames[concFieldInd]);
        OStringStream concFieldTitleStream;
        concFieldTitleStream << word(fieldFilePath.name()) << "_" << fieldFilePath.path().name();
        word concFieldTitle(concFieldTitleStream.str());
        ovFid<<"\t"<<concFieldTitle;
    }
    ovFid<<"\n";
    ovFid.close();

    qcFid.open(qcFileName.c_str());
    if(!qcFid.is_open())
      FatalErrorIn(args.executable())
	<< "Could not open file " << qcFileName << "for writing"
	<< exit(FatalError);
    
    qcFid<<"Quality assurance of calculations\n";
    qcFid.close();
  }
