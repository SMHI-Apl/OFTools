    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

/*
    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")       
    );
*/

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar PrT
    (
        transportProperties.lookup("PrT")
    );

    dimensionedScalar PrTt
    (
        transportProperties.lookup("PrTt")
    );
/*
    dimensionedScalar Cp
    (
        transportProperties.lookup("Cp")
    );

    dimensionedScalar Q
    (
        transportProperties.lookup("Q")
    );
*/
