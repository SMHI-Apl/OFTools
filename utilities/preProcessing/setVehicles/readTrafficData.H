
IOdictionary trafficDict
  (
    IOobject
    (
      "trafficDict",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

//Height of volume for emission release
//(
// readScalar(trafficData.lookup("vehicleHeight"))
//);

scalar vehicleArea
(
 readScalar(trafficDict.lookup("Acar"))
);
dimensionedScalar Acar("Acar",dimArea,vehicleArea);
//Acar=0;


wordList roadNames;
trafficDict.lookup("roadNames") >> roadNames;

//Reading road properties: speed, nr lanes in road, ADT
List< FixedList<scalar,3> > roadProperties;
trafficDict.lookup("roadProperties") >> roadProperties;


List< List<vector> >	centreLines;
trafficDict.lookup("centreLines") >> centreLines;

