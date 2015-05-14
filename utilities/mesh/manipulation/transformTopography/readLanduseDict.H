
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


scalar subtractedX(readScalar(landuseDict.lookup("subtractedX")));
scalar subtractedY(readScalar(landuseDict.lookup("subtractedY")));
