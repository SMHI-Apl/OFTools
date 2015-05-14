    Info << "\nReading environmentalProperties" << endl;

    IOdictionary environmentalProperties
    (
        IOobject
        (
            "environmentalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedVector g(environmentalProperties.lookup("g"));

    dimensionedScalar Omega(environmentalProperties.lookup("Omega"));

    dimensionedVector OmegaProjection(environmentalProperties.lookup("OmegaProjection"));

    dimensionedVector twoOmega(2*Omega*OmegaProjection);



