    Info<< "Reading lowerBoundValues\n" << endl;

    IOdictionary lowerBoundValues
    (
        IOobject
        (
            "lowerBoundValues",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar TMin
    (
        lowerBoundValues.lookup("TMin")
    );

    dimensionedScalar SMin
    (
        lowerBoundValues.lookup("SMin")
    );

