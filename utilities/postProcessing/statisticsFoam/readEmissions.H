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

 // air density to calculate emission as a massfraction instead of mass/s 
 dimensionedScalar densityAir(emissionDict.lookup("densityAir")); 

 // time-scale in seconds for chemical transformation (average transport time from source)
 scalar tau(readScalar(emissionDict.lookup("tau")));


 // Fraction of NOx directly emitted as NO2 (~10 %)
 scalar fd(readScalar(emissionDict.lookup("fd")));

