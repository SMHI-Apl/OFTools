  IOdictionary setForcesDict
  (
    IOobject
    (
      "setForcesDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );


Info<< "Reading set names from dictionary" << endl;
wordList setNames;
setForcesDict.lookup("setNames") >> setNames;
