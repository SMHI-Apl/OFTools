
   // Solve the Momentum equation

//volVectorField Urel;
    tmp<fvVectorMatrix> UEqn;
    if(traffic)
      {
	//Urel=U-Ucar;
	UEqn=
	  (
	   fvm::div(phi, U)
	   + turbulence->divDevReff(U)
	   ==
	   - 0.5*Cfcar*VAD*mag(U-Ucar)*(U-Ucar)
	   //	   -fvm::Sp(0.5*Cfcar*VAD*mag(Urel),Urel)
	   );

      }
    else 
      {
	UEqn=
	  (
	   fvm::div(phi, U)
	   + turbulence->divDevReff(U)
	   );
      }

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
