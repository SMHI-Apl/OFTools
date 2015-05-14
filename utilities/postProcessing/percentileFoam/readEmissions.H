IOdictionary emissionDict
  (
    IOobject
    (
      "emissionDict",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

 //air density to calculate emission as a massfraction instead of mass/s 
 dimensionedScalar densityAir(emissionDict.lookup("densityAir")); 
