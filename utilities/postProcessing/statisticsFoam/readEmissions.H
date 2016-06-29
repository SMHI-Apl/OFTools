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

 scalar tau=200;
 scalar fd=0.1;

if(NO2_conversion) {
    Info << "NO2-conversion activated" << endl;
    // time-scale in seconds for chemical transformation
    // (average transport time from source)
    tau = readScalar(emissionDict.lookup("tau"));
    Info << "Mixing time-scale is: " << tau << endl;
    
    // Fraction of NOx directly emitted as NO2 (~10 %)
    fd = readScalar(emissionDict.lookup("fd"));
    Info << "Fraction of direct emitted NO2 is: " << fd << endl;
}
