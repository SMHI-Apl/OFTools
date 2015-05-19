    Info<< "Creating landuse field" << endl;
    volScalarField landuse
    (
        IOobject
        (
            "landuse",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("landuse",dimensionSet (0,0,0,0,0,0,0),-1)
    );

    Info<< "Reading boundaryfield nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("nut",dimensionSet (0,2,-1,0,0,0,0),0)
    );

    Info<< "Creating field LAI\n" << endl;
    volScalarField LAD
    (
        IOobject
        (
            "LAD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("LAD",dimensionSet (0,-1,0,0,0,0,0),0)
    );

    // Info<< "Reading boundaryfield U\n" << endl;
    // volVectorField U
    // (
    //     IOobject
    //     (
    //         "U",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // Info << "Reading field epsilon\n" << endl;
    // volScalarField epsilon
    // (
    //     IOobject
    //     (
    //         "epsilon",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // Info << "Reading field k\n" << endl;
    // volScalarField k
    // (
    //     IOobject
    //     (
    //         "k",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );
