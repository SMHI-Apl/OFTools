    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
      - fvm::Sp(canopy*mag(U),U)
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
