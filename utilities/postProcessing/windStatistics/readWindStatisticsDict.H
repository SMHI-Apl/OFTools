Info<< "Reading wind statistics dictionary\n" << endl;

  IOdictionary windStatisticsDict
  (
    IOobject
    (
      "windStatisticsDict",
      runTime.system(),
      runTime,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  IOdictionary landuseDict
  (
    IOobject
    (
      "landuseDict",
      runTime.constant(),
      runTime,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

scalar subtractedX = readScalar(landuseDict.lookup("subtractedX"));
scalar subtractedY = readScalar(landuseDict.lookup("subtractedY"));

fileName mesanLine;
windStatisticsDict.lookup("mesanLine") >> mesanLine;
timeSeries mesanLineTS(mesanLine.c_str());
Info <<"Read mesanLineTS"<<endl;
int TInd=mesanLineTS.getColInd("T");
int cloudInd=mesanLineTS.getColInd("clouds");
if(TInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter T not found in timeseries header"<< exit(FatalError);

if(cloudInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter clouds not found in timeseries header"<< exit(FatalError);
mesanLineTS.removeNodata(TInd);
mesanLineTS.removeNodata(cloudInd);

fileName waspVanern;
windStatisticsDict.lookup("waspVanern") >> waspVanern;
timeSeries waspVanernTS(waspVanern.c_str());
Info <<"Read waspVanernTS"<<endl;
int waspVanernWspeedInd=waspVanernTS.getColInd("wspeed");
if(waspVanernWspeedInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter wspeed not found in timeseries waspVanern header"<< exit(FatalError);
waspVanernTS.removeNodata(waspVanernWspeedInd);

fileName mesanVanern;
windStatisticsDict.lookup("mesanVanern") >> mesanVanern;
timeSeries mesanVanernTS(mesanVanern.c_str());
Info <<"Read mesanVanernTS"<<endl;
int mesanVanernWdirInd=mesanVanernTS.getColInd("wdir");
if(mesanVanernWdirInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter wdir not found in mesanVanern timeseries header"<< exit(FatalError);
mesanVanernTS.removeNodata(mesanVanernWdirInd);
int mesanVanernWspeedInd=mesanVanernTS.getColInd("wspeed");
if(mesanVanernWspeedInd==-1)
    FatalErrorIn(args.executable())
      << "Parameter wspeed not found in mesanVanern timeseries header"<< exit(FatalError);
mesanVanernTS.removeNodata(mesanVanernWspeedInd);

fileName waspCoast;
windStatisticsDict.lookup("waspCoast") >> waspCoast;
timeSeries waspCoastTS(waspCoast.c_str());
Info <<"Read waspCoastTS"<<endl;
int waspCoastWspeed100Ind=waspCoastTS.getColInd("wspeed100");
int waspCoastWdirInd=waspCoastTS.getColInd("wdir");
if(waspCoastWspeed100Ind==-1)
    FatalErrorIn(args.executable())
      << "Parameter wspeed100 not found in waspCoast timeseries header"<< exit(FatalError);
waspCoastTS.removeNodata(waspCoastWspeed100Ind);

Info<<"Matching dates..."<<endl;
//waspVanernTS.matchDates(mesanVanernTS);
//waspVanernTS.matchDates(waspCoastTS);
//waspVanernTS.matchDates(mesanLineTS);

Info<<"Matching dates...one done"<<endl;
//mesanVanernTS.matchDates(waspVanernTS);
mesanVanernTS.matchDates(waspCoastTS);
mesanVanernTS.matchDates(mesanLineTS);

Info<<"Matching dates...two done"<<endl;
waspCoastTS.matchDates(mesanVanernTS);
//waspCoastTS.matchDates(waspVanernTS);
waspCoastTS.matchDates(mesanLineTS);

Info<<"Matching dates...three done"<<endl;
mesanLineTS.matchDates(mesanVanernTS);
mesanLineTS.matchDates(waspCoastTS);
//mesanLineTS.matchDates(waspVanernTS);
Info<<"Matching dates...all done"<<endl;

Info<<"Hours to be used for statistics: "<<mesanLineTS.nrows<<endl;

SortableList<label> wdirs(windStatisticsDict.lookup("wdirs"));
wdirs.sort();

vector refPoint;
windStatisticsDict.lookup("refPoint") >> refPoint;
scalar refHeight = readScalar(windStatisticsDict.lookup("rasterRefHeight"));
scalar refz0=readScalar(windStatisticsDict.lookup("refz0"));

fileName dirName;
windStatisticsDict.lookup("rasterArchiveDirectory") >> dirName;
 
fileName resultDir;
windStatisticsDict.lookup("resultDirectory") >> resultDir;
  
PtrList<Raster> rasters(wdirs.size());
List<scalar> refWinds(wdirs.size());


List<vector> points;
windStatisticsDict.lookup("points") >> points;


List<vector> probePointList;
windStatisticsDict.lookup("probePoints") >>probePointList;
List<word> probeNameList;
windStatisticsDict.lookup("probePointNames")>>probeNameList;
fileName probeFileName;
windStatisticsDict.lookup("probeFileName") >> probeFileName;


forAll(wdirs,wdiri)
{
  Info<<"Processing wdir: "<<wdirs[wdiri]<<endl;
  rasters.set(wdiri, new Raster );
  OStringStream rasterNameStream;
  fileName rasterName;
  rasterNameStream <<"wdir_"<<wdirs[wdiri]<<".asc";
  rasterName=dirName/rasterNameStream.str();
  
  //fileName rasterName=dirName+"/wdir_"+word(wdirs[wdiri])+".asc";
  // rasterName <<dirName<<"/wdir_"<<wdirs[wdiri]<<".asc";
    //nameStr.str();
  if (!rasters[wdiri].read(rasterName.c_str()))
    FatalErrorIn(args.executable())
      << "Could not read raster file: "<<rasterName<< exit(FatalError);
  else
    Info<<"Read raster, max: "<<rasters[wdiri].max()<<"m/s, min: "<<rasters[wdiri].min()<<" m/s"<<endl;

  OStringStream refWindSS;
  refWindSS << "refWind."<<wdirs[wdiri];
  //word refWindLabel;
  //refWindLabel = refWindSS.str();
  scalar refWind=readScalar(windStatisticsDict.lookup(refWindSS.str()));
  Info<<"Ref wind on height"<<refPoint[2]<<": "<<refWind<<" m/s"<<endl;
  //  scalar ustar=0.41*refWind/(Foam::log((refPoint[2]+z0)/z0));
  refWind=refWind*Foam::log((refHeight+refz0)/refz0)/Foam::log((refPoint[2]+refz0)/refz0);
  refWinds[wdiri]=refWind;
  Info<<"Ref wind on raster height: "<<refWind<<" m/s"<<endl;
}


  
