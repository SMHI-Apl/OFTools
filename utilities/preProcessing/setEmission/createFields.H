
// List of source term fields
 PtrList<volScalarField> eST(emissionGroups.size());

 word fieldPrefix = "eST_";
 word fieldTitle;
 
 for(label i=0;i<emissionGroups.size(); i++)
      {
	fieldTitle = fieldPrefix & emissionGroups[i];
	eST.set
              (
	          i,
                  new volScalarField
                  (
                      IOobject
                      (
                          fieldTitle,
                          runTime.timeName(),
                          mesh,
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE
                      ),
                      mesh,
		      dimensionedScalar("eST_init",dimensionSet(1,-3,-1,0,0,0,0),0)
                  )
              );
      }

volScalarField eSTij
(
 IOobject
 (
  "eSTij",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 eST[0]
 );
