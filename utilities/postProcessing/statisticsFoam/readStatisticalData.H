Info<< "Reading statistical data\n" << endl;

  IOdictionary statisticsDict
  (
    IOobject
    (
      "statisticsDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  Switch write3DFields(statisticsDict.lookup("write3DFields"));
  scalar hourlyPercentile(readScalar(statisticsDict.lookup("hourlyPercentile")));
  scalar dailyAvgPercentile(readScalar(statisticsDict.lookup("dailyAveragePercentile")));
  scalar dailyMaxPercentile(readScalar(statisticsDict.lookup("dailyMaxPercentile")));
  
  scalar nValidsForDaily(readScalar(statisticsDict.lookup("nValidsForDaily")));
  
  scalar maxConc(readScalar(statisticsDict.lookup("maxAllowedConc")));
  Info<<"Max allowed concentration is: "<<maxConc<<" ug/m3"<<endl;
  fileName metFile;
  statisticsDict.lookup("metFile") >> metFile;
  
  fileName emisFile;
  statisticsDict.lookup("emisFile") >> emisFile;

  fileName probeDir;
  statisticsDict.lookup("probeDir") >> probeDir;

  Switch background(statisticsDict.lookup("backgroundSwitch"));
  Switch NO2_conversion(statisticsDict.lookup("NO2_conversionSwitch"));
  
  fileName bgFile;
  statisticsDict.lookup("backgroundFile") >> bgFile;

  fileName concFileDir;
  statisticsDict.lookup("concFileDir") >> concFileDir;
  List<fileName> concFiles;
  statisticsDict.lookup("concFileList")>>concFiles;

  std::vector<std::string> concFileNames;
  for(int i=0; i<int(concFiles.size());i++) {
      OStringStream concFileNameStream;
      concFileNameStream << word(concFileDir) << word(concFiles[i]);
      concFileNames.push_back(std::string(concFileNameStream.str()));
  }

  List<vector> probeList;
  statisticsDict.lookup("probeLocations")>> probeList;

  wordList probeNames;
  statisticsDict.lookup("probeNames")>> probeNames;

  if(probeList.size() != probeNames.size())
    FatalErrorIn(args.executable())
        << "Number of probes does not match number of probe names given in dictionary: statisticsDict"
        << exit(FatalError);

  Info<<"Statistical data has been read\n" << endl;

  if(concFiles.size()==0)
    FatalErrorIn(args.executable())
        << "No concentration files specified in dictionary: statisticsDict"
        << exit(FatalError);
  else
    Info<<"Number of concentration files given in statisticsDict is:"
        <<int(concFileNames.size())<<endl;

  Info<< "Inititalizing percentile calculation!"<<endl;
  Info<< "Percentile of hourly concentrations: "<<hourlyPercentile<<endl;
  Info<< "Percentile of daily average concentrations: "<<dailyAvgPercentile<<endl;
  Info<< "Percentile of daily max concentrations: "<<dailyMaxPercentile<<endl;
  Info<< "Daily values are based on at least "<<nValidsForDaily<<" hourly values"<<endl;

  Info<< "Meteorological data from: "<<metFile<<endl;
  Info<< "Emissions from: "<<emisFile<<endl;
  if(background)
      Info<< "Background concentrations from: "<<bgFile<<endl;
  else
      Info<< "No background concentrations requested"<<endl;
