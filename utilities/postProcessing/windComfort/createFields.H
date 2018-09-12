

//Used later to avoid devide by zero in analysis
volScalarField Umin
 (
  IOobject
  (
   "Umin",
   runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::NO_WRITE
   ),
  mesh,
  dimensionedScalar("Umin",dimensionSet (0,1,-1,0,0,0,0),0.01)
  );

 volScalarField UekvMedian
 (
  IOobject
  (
   "UekvMedian",
   runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::AUTO_WRITE
   ),
  mesh,
  dimensionedScalar("Uekv",dimensionSet (0,1,-1,0,0,0,0),0)
  );


 volScalarField UekvProc
 (
  IOobject
  (
   "UekvProc",
   runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::AUTO_WRITE
   ),
  mesh,
  dimensionedScalar("proc",dimensionSet (0,0,0,0,0,0,0),0)
  );


