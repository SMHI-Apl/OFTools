    


// Solving pollutant specie transport equation
// Solving for n-1 species
for(label i = 0; i < Y.size(); i++)
  {      
    volScalarField& Yi = Y[i];
    volScalarField& eSTi = eST[i];
    Info<<"hej1"<<endl;
    //solve
    //  (
    //   fvm::div(phi, Yi, divScheme)
    //   - fvm::laplacian((turbulence->nut())/SC_T, Yi, laplacianScheme)
    //   - eSTi/densityAir,
    //   // Solver
    //   mesh.solver("Yi")
    //   );

    solve
      (
       fvm::div(phi, Yi, divScheme)
       - fvm::laplacian((turbulence->nut())/SC_T, Yi, laplacianScheme)
       ==
       eSTi/densityAir
       );
    Info<<"hej2"<<endl;
    // mesh.solver("Yi")    
    Yi.max(0.0);
    //Info<< "Done specie: "<< i << endl;        
  }


