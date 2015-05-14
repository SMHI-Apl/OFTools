    Info << "\nReading coriolisProperties" << endl;

    IOdictionary coriolisProperties
    (
        IOobject
        (
            "coriolisProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Omega(coriolisProperties.lookup("Omega"));

    dimensionedVector OmegaProjection(coriolisProperties.lookup("OmegaProjection"));

    dimensionedVector twoOmega(2.0*Omega*OmegaProjection);
