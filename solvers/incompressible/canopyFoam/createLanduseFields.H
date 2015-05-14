
IOdictionary landuseDict
(
 IOobject
 (
  "landuseDict",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

Info<< "Reading landuse types from dictionary" << endl; 
tensorField landuseList(landuseDict.lookup("landuseList"));

wordList sourcePatches(landuseDict.lookup("sourcePatches"));

volScalarField landuse
(
 IOobject
 (
  "landuse",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField LAD
(
 IOobject
 (
  "LAD",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

volScalarField canopy
(
 IOobject
 (
  "canopy",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("canopy",dimensionSet(0,-1,0,0,0,0,0),0)
 );

HashTable<tensor,label> codeTable;
forAll(landuseList,codei)
{
  tensor code = landuseList[codei];
  codeTable.insert(label(code.xx()),code);
}

forAll(landuse.internalField(),celli)
{

  tensor code = codeTable[label(landuse.internalField()[celli])];
  scalar Cd = code.xy();
  scalar fraction = code.yx(); 
  canopy.internalField()[celli]=Cd*fraction*LAD.internalField()[celli];
}

forAll(sourcePatches,sp)
{
  label patchI=mesh.boundaryMesh().findPatchID(sourcePatches[sp]);
  
  forAll(canopy.boundaryField()[patchI],facei)
    {
      tensor code = codeTable[label(landuse.boundaryField()[patchI][facei])];
      canopy.boundaryField()[patchI][facei] = code.yy();
    }
}	  

canopy.write();
Info<<"Wrote canopy (fix to make field accesible to turbulence model"<<endl;
