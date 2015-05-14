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

    dimensionedScalar T0
    (
        lowerBoundValues.lookup("T0")
    );

    dimensionedScalar S0
    (
        lowerBoundValues.lookup("S0")
    );

