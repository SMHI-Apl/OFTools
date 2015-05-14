{
  //    volScalarField kappatS("kappatS",turbulence->nut()/PrtS);
    volScalarField kappatS=turbulence->nut()/PrtS;
    //    kappatS.correctBoundaryConditions(); 

    volScalarField kappaEffS("kappaEffS", turbulence->nu()/PrS + kappatS);
    kappatS.clear();

    fvScalarMatrix SEqn
    (
        fvm::ddt(S)
      + fvm::div(phi, S)
      - fvm::laplacian(kappaEffS, S)
    );

    kappaEffS.clear();

    SEqn.relax();
    SEqn.solve();

    //    rhok = 1.0 - beta*(T - TRef);
	
    volScalarField TC=T-T0C;//;273.15;

    rhok = (a0 + (a1 + (a2 + (a3 + (a4 + a5*TC)*TC)*TC)*TC)*TC
      + (b0 + (b1 + (b2 + (b3 + b4*TC)*TC)*TC)*TC)*S
	    + (c0 + (c1 + c2*TC)*TC)*pow(S,3/2) + d0*S*S)/a0;

    TC.clear();
}
