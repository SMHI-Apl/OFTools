Info<< "Reading wind comfort dictionary\n" << endl;

  IOdictionary windComfortDict
  (
    IOobject
    (
      "windComfortDict",
      runTime.system(),
      runTime,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

scalar criteriaUlim = readScalar(windComfortDict.lookup("criteriaUlim"));
scalar Iref = readScalar(windComfortDict.lookup("Iref"));
scalar alpha = readScalar(windComfortDict.lookup("alpha"));
List<vector> refPoints(windComfortDict.lookup("refPoints"));
scalar metTSNodata = readScalar(windComfortDict.lookup("metTSNodata"));

fileName metFileName;

windComfortDict.lookup("metFileName") >> metFileName;
timeSeries metTS(metFileName.c_str());
Info <<"Successfully read meteorology time-series"<<endl;
scalar detectionLimit =readScalar(windComfortDict.lookup("detectionLimit"));
scalar wspeed = readScalar(windComfortDict.lookup("wspeed")); //To open archived fields

//fileName dirName;
//windComfortDict.lookup("flowArchiveDirectory") >> dirName;
 
List<label> wdirsTmp(windComfortDict.lookup("wdirs"));
SortableList<label> wdirs(wdirsTmp);
//wdirs.sort();
if(wdirs.size()<2)
    FatalErrorIn(args.executable())
      << "At least two wind directions needed"<< exit(FatalError);

int wspeedInd=metTS.getColInd("wspeed");
int wdirInd=metTS.getColInd("wdir");
if(wspeedInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter 'wspeed' not found in timeseries header"<< exit(FatalError);

if(wdirInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter 'wdir' not found in timeseries header"<< exit(FatalError);

//All hours with unknown windspeed are deleted
//Calm weather should be marked with something else than nodata
metTS.removeNodata(double(metTSNodata),wspeedInd);

Info<<"Hours to be used for statistics: "<<metTS.nrows<<endl;
Info<<endl;

  






