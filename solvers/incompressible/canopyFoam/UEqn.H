    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        - fvm::Sp(canopy*mag(U),U) + fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
