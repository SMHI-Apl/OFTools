  IOdictionary createFaceSetsDict
  (
    IOobject
    (
      "createFaceSetsDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );


Info<< "Reading patch names from dictionary" << endl;
wordList patchNames;
createFaceSetsDict.lookup("patchNames") >> patchNames;

Info<< "Reading limit values from dictionary" << endl;
scalarList xlims,ylims,zlims;
createFaceSetsDict.lookup("xlims")>>xlims;
createFaceSetsDict.lookup("ylims")>>ylims;
createFaceSetsDict.lookup("zlims")>>zlims;
Info<<"xlims: "<<xlims<<endl;
Info<<"ylims: "<<ylims<<endl;
Info<<"zlims: "<<zlims<<endl;

vector normal;
createFaceSetsDict.lookup("normal") >> normal;
scalar tol(readScalar(createFaceSetsDict.lookup("tol")));
word prefix;
createFaceSetsDict.lookup("prefix") >> prefix;
