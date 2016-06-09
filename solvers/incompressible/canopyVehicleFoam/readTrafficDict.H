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


Switch traffic(trafficDict.lookup("traffic"));

scalar Cfcar(readScalar(trafficDict.lookup("Cfcar")));
