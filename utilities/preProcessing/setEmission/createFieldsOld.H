
// List of source term fields
 PtrList<volScalarField> eST(emissionGroups.size());

 word fieldPrefix = "eST_";
 word fieldTitle;
 
 for(label i=0;i<emissionGroups.size(); i++)
      {
	fieldTitle = fieldPrefix & emissionGroups[i];
	  IOobject header
          (
              fieldTitle,
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ
          );

          // Checking if the field already exists and can be read
          if (header.headerOk())
          {
              eST.set
              (
	          i,
                  new volScalarField
                  (
                      IOobject
                      (
                          fieldTitle,
                          mesh.time().timeName(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::AUTO_WRITE
                      ),
                      mesh
                  )
              );
          }
          else
          {
              volScalarField eSTdefault
              (
                  IOobject
                  (
                      "eSTdefault",
                      mesh.time().timeName(),
                      mesh,
                      IOobject::MUST_READ,
                      IOobject::NO_WRITE
                  ),
                  mesh
              );

              eST.set
              (
	          i,
                  new volScalarField
                  (
                      IOobject
                      (
                          fieldTitle,
                          mesh.time().timeName(),
                          mesh,
                          IOobject::NO_READ,
                          IOobject::AUTO_WRITE
                      ),
                      eSTdefault
                  )
              );
          }
      }
volScalarField eSTij
(
 IOobject
 (
  "eSTij",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 eST[0]
 );
