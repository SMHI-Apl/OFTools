    singlePhaseTransportModel laminarTransport(U, phi);

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Laminar Prandtl number for salt transport equation
    dimensionedScalar PrS(laminarTransport.lookup("PrS"));

    // Turbulent Prandtl number for salt transport equation
    dimensionedScalar PrtS(laminarTransport.lookup("PrtS"));

    // T0C
    dimensionedScalar T0C(laminarTransport.lookup("T0C"));

