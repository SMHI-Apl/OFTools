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


wordList emissionGroupNames;
trafficDict.lookup("emissionGroupNames") >> emissionGroupNames;

wordList roadNames;
trafficDict.lookup("roadNames") >> roadNames;

List< FixedList<scalar,3> > roads;
trafficDict.lookup("roadClassification") >> roads;

// Table of species (in this case roads)
speciesTable emissionGroups(emissionGroupNames);


