    Info<< "Reading thermophysical properties\n" << endl;

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field S\n" << endl;
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    #include "myReadTransportProperties.H"
    #include "myReadDensityCoeffs.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
	((a0 + a1*(T-T0C) + a2*pow((T-T0C),2) + a3*pow((T-T0C),3) 
	+ a4*pow((T-T0C),4) + a5*pow((T-T0C),5))
	+ (b0 + b1*(T-T0C) + b2*pow((T-T0C),2) + b3*pow((T-T0C),3) + b4*pow((T-T0C),4))*S 
	+ (c0 + c1*(T-T0C) + c2*pow((T-T0C),2))*pow(S,3/2)
	+ d0*pow(S,2))/a0
    );

    // kinematic turbulent thermal thermal conductivity m2/s
    Info<< "Reading field kappat\n" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

