    Info<< "Reading emissionDict\n" << endl;

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

    // List of the species names(in this case road sources)
    wordList speciesNames(trafficDict.lookup("emissionGroupNames"));

    // Table of species
    speciesTable species(speciesNames);
    
    //air density to calculate emission as a massfraction instead of mass/s 
    dimensionedScalar densityAir(emissionDict.lookup("densityAir")); 

    // Turbulent Schmidt number
    dimensionedScalar SC_T
    (
        emissionDict.lookup("SC_T")
    );

    // Species fractions
    PtrList<volScalarField> Y(species.size());
    word speciesPrefix ="spec_";
    

   //Emission source terms
    PtrList<volScalarField> eST(species.size());    
    word sourceTermPrefix ="eST_";

    word fieldTitle;


   //********Creating species fractions scalar fields**********
   forAll(species, i)
    {
      fieldTitle=speciesPrefix & species[i];
      IOobject header
	(
	 fieldTitle,
	 runTime.timeName(),
	 mesh,
	 IOobject::NO_READ
	 );
           
      // Checking if the field already exists and can be read
      if (header.headerOk())
	{
	  Y.set
	    (
	     i,
	     new volScalarField
	     (
	      IOobject
	      (
	       fieldTitle,
	       runTime.timeName(),
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
	  volScalarField spec_default
	    (
	     IOobject
	     (
	      "spec_default",
	      runTime.timeName(),
	      mesh,
	      IOobject::MUST_READ,
	      IOobject::NO_WRITE
	      ),
	     mesh
	     );
	  
	  Y.set
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
		spec_default
		)
	       );
	}	  
    }


    //*******creating source term fields**********
    forAll(species, i)
     {
       fieldTitle=sourceTermPrefix & species[i];
       IOobject eST_header
	 (
	  fieldTitle,
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ
          );
       
       // Checking if the field already exists and can be read
       if (eST_header.headerOk())
	 {
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
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
		),
	       mesh
	       )
              );
	 }
       else
	 {
	   Info<<"Warning: No source term set for species "<<species[i]<<" using default zero source term"<<endl;
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
     }

				       
				       
   Info << "Reading field p\n" << endl;
   volScalarField p
    (
     IOobject
     (
      "p",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
      ),
     mesh
    );

   Info << "Reading field U\n" << endl;
   volVectorField U
    (
     IOobject
     (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
      ),
     mesh
     );

#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
     incompressible::RASModel::New(U, phi, laminarTransport)
    );

