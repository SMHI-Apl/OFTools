Info<< "Reading windDict\n" << endl;

  IOdictionary windDict
  (
    IOobject
    (
      "windDict",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  scalar direction(readScalar(windDict.lookup("windDirection")));
  scalar z0(readScalar(windDict.lookup("z0")));
  scalar Umag(readScalar( windDict.lookup("U10")));
  scalar top(readScalar( windDict.lookup("top")));
  scalar d(readScalar( windDict.lookup("displacement")));
  //scalar initK(readScalar( windDict.lookup("initK")));
  //scalar initEps(readScalar( windDict.lookup("initEps")));

Info<< "windDict has been read\n" << endl;
