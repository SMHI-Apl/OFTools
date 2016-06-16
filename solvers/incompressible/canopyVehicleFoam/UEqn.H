    // Momentum predictor

    tmp<fvVectorMatrix> UEqn;
if(traffic)
{
    UEqn=
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        - 0.5*Cfcar*VAD*mag(U-Ucar)*(U-Ucar)
        - fvm::Sp(canopy*mag(U),U) + fvOptions(U)
    );
}
else
{
    UEqn=
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        - fvm::Sp(canopy*mag(U),U) + fvOptions(U)
    );
}

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
