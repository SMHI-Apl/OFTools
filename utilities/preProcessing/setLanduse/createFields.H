 Info<< "Reading landuse field" << endl;
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
