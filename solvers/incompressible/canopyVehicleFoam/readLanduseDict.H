
  IOdictionary landuseDict
  (
    IOobject
    (
      "landuseDict",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

Switch canopy(landuseDict.lookup("canopy"));


Info<< "Reading source face names from dictionary" << endl;
wordList sourcePatches;
landuseDict.lookup("sourcePatches") >> sourcePatches;

Info<< "Reading landuse types from dictionary" << endl;
tensorField landuseList(landuseDict.lookup("landuseList")); 








