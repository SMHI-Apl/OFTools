    Info<< "Reading probeLocations\n" << endl;

    IOdictionary probeLocations
    (
        IOobject
        (
            "probeLocations",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Constructing probes\n" << endl;

    PtrList< Probe<volScalarField> >  probes(species.size());
    for(label i = 0; i < Y.size(); i++)
    {      
      //     volScalarField& Yi = Y[i];
      probes.hook
	(
	 new Probe<volScalarField>(probeLocations, Y[i])
	 );
    }

