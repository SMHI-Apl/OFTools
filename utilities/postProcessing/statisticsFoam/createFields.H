
int indHourlyPerc=0,indDailyAvgPerc=0, indDailyMaxPerc=0;

label nStatFields=1;

if(hourlyPercentile>0 and hourlyPercentile <100)
  {
    nStatFields++;
    indHourlyPerc=nStatFields-1;
  }

if(dailyAvgPercentile>0 and dailyAvgPercentile <100)
{
  nStatFields++;
  indDailyAvgPerc=nStatFields-1;;
}

if(dailyMaxPercentile>0 and dailyMaxPercentile <100)
{
  nStatFields++;
  indDailyMaxPerc=nStatFields-1;
}

// Number of concentration fields
label nConcFields = ca.fileNames.size();
PtrList<volScalarField> concFields(nConcFields);
if (args.parRunControl().parRun()) {
    Info<<"Parallel run, reading all concentration fields"<<endl;
    for(label concFieldInd=0; concFieldInd < nConcFields; concFieldInd++) {
        fileName fieldFilePath(ca.fileNames[concFieldInd]);
        fileName archiveDirPath = fieldFilePath.path();
        fileName fieldFileName = fieldFilePath.name();
        OStringStream concFieldTitleStream;
        concFieldTitleStream << word(fieldFileName) << "_" << archiveDirPath.name();
        word concFieldTitle(concFieldTitleStream.str());
        Info<<"Reading field " << concFieldTitle.c_str()<<endl;
        
        concFields.set(
            concFieldInd,
            new volScalarField
            (
                IOobject
                (
                    concFieldTitle,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }
}


PtrList<volScalarField> statFields(nStatFields);
label statFieldInd=0;
Info<< "Creating field for average concentrations"<<endl;
statFields.set(
	       statFieldInd,
	       new volScalarField
	       (
		IOobject
		(
		 "hourlyAvg",
		 runTime.timeName(),
		 mesh,
		 IOobject::NO_READ,
		 IOobject::AUTO_WRITE
		 ),
		mesh,
		dimensionedScalar("init_stat",dimensionSet(0,0,0,0,0,0,0),0)
		)
	       );
statFieldInd++;

if(indHourlyPerc)
{
  Info<< "Creating field for percentile of hourly concentrations"<<endl;
  std::stringstream titleStream;
  titleStream<<string("hourlyPerc_")<<hourlyPercentile;
  std::string titleString;
  titleStream>>titleString;

  word fieldTitle(titleString.c_str());
 
  statFields.set
    (
     statFieldInd,
     new volScalarField
     (
      IOobject
      (
       fieldTitle,
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
       ),
      mesh,
      dimensionedScalar("init_stat",dimensionSet(0,0,0,0,0,0,0),0)
      )
     );
  statFieldInd++;
}


if(indDailyAvgPerc)
{
 Info<< "Creating field for percentile of daily average concentrations"<<endl;
 std::stringstream titleStream;
 titleStream<<string("dailyAvgPerc_")<<dailyAvgPercentile;
 std::string titleString;
 titleStream>>titleString;
 
 word fieldTitle(titleString.c_str());
 
 statFields.set
   (
    statFieldInd,
    new  volScalarField
    (
     IOobject
     (
      fieldTitle,
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("init_stat",dimensionSet(0,0,0,0,0,0,0),0)
     )
    );
 
 statFieldInd++;
}

if(indDailyMaxPerc)
{
  Info<< "Creating field for percentile of daily max concentrations"<<endl;
  std::stringstream titleStream;
 titleStream<<string("dailyMaxPerc_")<<dailyMaxPercentile;
 std::string titleString;
 titleStream>>titleString;
 
 word fieldTitle(titleString.c_str());  

 statFields.set
   (
    statFieldInd,
    new volScalarField
    (
     IOobject
     (
      fieldTitle,
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar("init_stat",dimensionSet(0,0,0,0,0,0,0),0)
     )
    );
 statFieldInd++;
}

Info<<"All requested fields created!"<<endl;
 
