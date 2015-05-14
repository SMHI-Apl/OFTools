
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

scalarList patchLanduse;
landuseDict.lookup("patchLanduse") >> patchLanduse;
scalarList heightDistribution;
landuseDict.lookup("heightDistribution")>>heightDistribution;

Info<< "Reading landuse types from dictionary" << endl;
tensorField landuseList(landuseDict.lookup("landuseList")); 

scalar subtractedX(readScalar(landuseDict.lookup("subtractedX")));
scalar subtractedY(readScalar(landuseDict.lookup("subtractedY")));

Info<< "Mesh translation offsets are read from dictionary"<<endl;

Info<< "X offset is: "<<-1*subtractedX<<endl;
Info<< "Y offset is: "<<-1*subtractedY<<endl;
