 // List of source term fields
 PtrList<volScalarField> Udirs(wdirs.size());
 PtrList<volScalarField> kdirs(wdirs.size());


 Info<< "Reading flow fields"<<endl;
 forAll(wdirs,i)
     {
       
	OStringStream UfileNameStream;
	OStringStream kfileNameStream;

	UfileNameStream << "U_wspeed_" << fixed << setprecision(1);
	UfileNameStream << float(wspeed) << "_wdir_" << float(wdirs[i]);
	word UfieldTitle(UfileNameStream.str());
	Info <<"Reading field "<<UfieldTitle<<endl;
	volVectorField U
	(
	 IOobject
	 (
	  UfieldTitle,
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh
	);

	volScalarField Umag = mag(U);

	forAll(Umag, celli)
	  {
	    if(Umag[celli]<0.01)
	      Umag[celli]=0.01;
	  }
	Udirs.set(i, Umag);

	kfileNameStream << "k_wspeed_" << fixed << setprecision(1);
	kfileNameStream << float(wspeed) << "_wdir_" << float(wdirs[i]);
	word kfieldTitle(kfileNameStream.str());
	Info <<"Reading field "<<kfieldTitle<<endl;
	volScalarField k
	(
	 IOobject
	 (
	  kfieldTitle,
	  runTime.timeName(),
	  mesh,
	  IOobject::MUST_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh
	);

	kdirs.set(i, k);
      }

Info << "Size of Udirs is: " << Udirs.size() << endl;

