
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

Info<< "Reading source face names from dictionary" << endl;
wordList sourcePatches;
landuseDict.lookup("sourcePatches") >> sourcePatches;


