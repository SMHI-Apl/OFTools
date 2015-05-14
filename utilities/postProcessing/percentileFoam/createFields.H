
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

PtrList<volScalarField> statFields(nStatFields);

Info<< "Creating field for average concentrations"<<endl;
label fieldInd=0;
statFields.set(
	       fieldInd,
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
fieldInd++;

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
     fieldInd,
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
  fieldInd++;
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
    fieldInd,
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
 
 fieldInd++;
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
    fieldInd,
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
 fieldInd++;
}

Info<<"All requested fields created!"<<endl;
 
