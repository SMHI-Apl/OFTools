    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
        ==
        flowDirection*gradP
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
