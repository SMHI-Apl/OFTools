    Info<< "\nReading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );


    //  Read centerline velocity for channel simulations
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );

    dimensionedScalar magUbar = mag(Ubar);
    vector flowDirection = (Ubar/magUbar).value();
