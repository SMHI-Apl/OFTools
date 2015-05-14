

    vector Ucar_init;
    Ucar_init[0]=0;
    Ucar_init[1]=0;
    Ucar_init[2]=0;

    // Car average speed
    volVectorField Ucar
    (
        IOobject
        (
            "Ucar",
            word("0"),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedVector(word("Ucar"),dimensionSet(0, 1, -1, 0, 0, 0, 0),Ucar_init)
    );


    //Vehicle Area Density
    volScalarField VAD
    (
        IOobject
        (
            "VAD",
            word("0"),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(word("VAD"),dimensionSet(0, -1, 0, 0, 0, 0, 0),0)
    );


    volScalarField fluidFrac
    (
        IOobject
        (
            "fluidFrac",
            word("0"),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(word("fluidFrac"),dimensionSet(0, 0, 0, 0, 0, 0, 0),1.0)
    );
