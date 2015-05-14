    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    Info << "Reading field U\n" << endl;
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

#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

/*
    Info<< "Reading/calculating field TC\n" << endl;
    volScalarField TC
    (
        IOobject
        (
            "TC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	T - T0C
    );
*/
    Info<< "Reading/calculating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
	(a0 + a1*(T-T0C) + a2*pow((T-T0C),2) + a3*pow((T-T0C),3) 
        + a4*pow((T-T0C),4) + a5*pow((T-T0C),5))
	+ (b0 + b1*(T-T0C) + b2*pow((T-T0C),2) + b3*pow((T-T0C),3) + b4*pow((T-T0C),4))*S 
	+ (c0 + c1*(T-T0C) + c2*pow((T-T0C),2))*pow(S,3/2)
	+ d0*pow(S,2)
    );

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
