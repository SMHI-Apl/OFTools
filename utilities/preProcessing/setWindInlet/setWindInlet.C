/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    setWindInlet

Description
    Set the inlet velocities for wind calculations.

\*---------------------------------------------------------------------------*/

//include "database.H"

#include "fvCFD.H"
#include "OSspecific.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readWindDict.H"

//   include "createDatabase.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  const scalar pi=3.1415926536;
  scalar dirRad;
  vector U10(0, 0, 0);
  int number_of_patches=1; 

  //converts to radians and unit circle
  dirRad=pi/180*(90.-direction);
  U10.x()=-Umag*Foam::cos(dirRad);
  U10.y()=-Umag*Foam::sin(dirRad);

  if(U10.x()<1e-6 && U10.x()>-1e-6)
    U10.x()=0.0;
  if(U10.y()<1e-6 && U10.y()>-1e-6)
    U10.y()=0.0;
  
  Info<< "Wind mag. is " << Umag << " ,Wind dir. is "<< direction <<" degrees" << endl;
  Info<< "X comp. of U at 10 m is " << U10.x() << endl;
  Info<< "Y comp. of U at 10 m is " << U10.y() << endl;

  //Setting initial internal values
  // k.internalField()=initK;
  //epsilon.internalField()=initEps;
  U.internalField()=U10;


  // Declaring patch index labels
  label patchI;
  label patchII;

  // Finding indexes for inlet patches corresponding to given wind direction
  if (direction==360 || direction==0){
    patchI = mesh.boundaryMesh().findPatchID("north");
    Info<< "Found patch north at index " << patchI << endl;}
  else if (direction == 90) {
    patchI = mesh.boundaryMesh().findPatchID("east");
    Info<< "Found patch east at index " << patchI << endl;}
  else if (direction == 180){
    patchI = mesh.boundaryMesh().findPatchID("south");
    Info<< "Found patch south at index " << patchI << endl;}  
  else if (direction == 270){
    patchI = mesh.boundaryMesh().findPatchID("west");
    Info<< "Found patch west at index " << patchI << endl;}
  else 
    {
      number_of_patches=2;
      if(direction>0 && direction <90)
	{     
	  patchI=mesh.boundaryMesh().findPatchID("north");
	  patchII=mesh.boundaryMesh().findPatchID("east");
	  Info<< "Found patches north and east at index " << patchI <<" and" << patchII<< endl; 
	}
      else if(direction>90 && direction <180)
	{     
	  patchI=mesh.boundaryMesh().findPatchID("east");
	  patchII=mesh.boundaryMesh().findPatchID("south");
	  Info<< "Found patches east and south at index " << patchI <<" and" << patchII<< endl; 
	}
      else if(direction>180 && direction <270)
	{     
	  patchI=mesh.boundaryMesh().findPatchID("south");
	  patchII=mesh.boundaryMesh().findPatchID("west");
	  Info<< "Found patches south and west at index " << patchI <<" and" << patchII<< endl; 
	}
      else
	{     
	  patchI=mesh.boundaryMesh().findPatchID("west");
	  patchII=mesh.boundaryMesh().findPatchID("north");
	  Info<< "Found patches west and north at index " << patchI <<" and" << patchII<< endl; 
	}
      if (patchII == -1)
	{
	  Info<< "Cannot find inlet patch\n" << endl;
	  FatalErrorIn(args.executable())
	    << "Cannot find inlet patch" << exit(FatalError);
	}
    }
  
  if (patchI == -1)
    {
      Info<< "Cannot find inlet patch\n" << endl;
      FatalErrorIn(args.executable())
	<< "Cannot find inlet patch" << exit(FatalError);
    }   
  
  
  // Declaring and initializing references to boundary patches
  const polyPatch& pp = mesh.boundaryMesh()[patchI];
 

  // Declaring and initializing references to datafields
  fixedValueFvPatchVectorField& inletVelocity =
    refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchI]);  
  fixedValueFvPatchScalarField& inletTurb = 
    refCast<fixedValueFvPatchScalarField>(k.boundaryField()[patchI]);
  fixedValueFvPatchScalarField& inletEps = 
    refCast<fixedValueFvPatchScalarField>(epsilon.boundaryField()[patchI]);

  scalar kappa=0.41;
  scalar Cmy=0.09;
  scalar ustar=kappa*Umag/(Foam::log((10+z0)/z0));  
      
  // Looping over all faces of first patch and assigning boundary values
  forAll(inletVelocity, patchFaceI)
    {
      const point& faceCentre =pp.faceCentres()[patchFaceI];      
      scalar z = faceCentre.z();
    
      if ( (z-d)>0 && z<top )
	inletVelocity[patchFaceI] =  U10*Foam::log((z-d+z0)/z0)/Foam::log((10.+z0)/z0);
      else if (z-d<0)
	inletVelocity[patchFaceI] = U10*0.0;
      else
	inletVelocity[patchFaceI] =  U10*Foam::log((top-d+z0)/z0)/Foam::log((10.+z0)/z0);


      inletTurb[patchFaceI] = Foam::pow(ustar,2)/Foam::sqrt(Cmy);

      if(z-d+z0>0)
	inletEps[patchFaceI] = Foam::pow(ustar,3)/(kappa*(z-d+z0));      
      else
	inletEps[patchFaceI] = Foam::pow(ustar,3)/(kappa*(z0));
    }
   
  if(number_of_patches == 2)
    {
      const polyPatch& pp2= mesh.boundaryMesh()[patchII];
      // Declaring and initializing references to datafields
      fixedValueFvPatchVectorField& inletVelocity2 =
	refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchII]);  
      fixedValueFvPatchScalarField& inletTurb2 = 
	refCast<fixedValueFvPatchScalarField>(k.boundaryField()[patchII]);
      fixedValueFvPatchScalarField& inletEps2 = 
	refCast<fixedValueFvPatchScalarField>(epsilon.boundaryField()[patchII]);
      
      // Looping over all faces of first patch and assigning boundary values
      forAll(inletVelocity2, patchFaceI)
	{
	  const point& faceCentre =pp2.faceCentres()[patchFaceI];      
	  scalar z = faceCentre.z();
	  
	  if ( (z-d)>0 && z<top )
	    inletVelocity2[patchFaceI] =  U10*Foam::log((z-d+z0)/z0)/Foam::log((10.+z0)/z0);
	  else if (z-d<0)
	    inletVelocity2[patchFaceI] = U10*0.0;
	  else
	    inletVelocity2[patchFaceI] =  U10*Foam::log((top-d+z0)/z0)/Foam::log((10.+z0)/z0);
	  

	  inletTurb2[patchFaceI] = Foam::pow(ustar,2)/Foam::sqrt(Cmy);
	  if(z-d+z0>0)
	    inletEps2[patchFaceI] = Foam::pow(ustar,3)/(kappa*(z-d+z0));      
	  else
	    inletEps2[patchFaceI] = Foam::pow(ustar,3)/(kappa*(z0)); 
	}
    }
  scalar z;
  vector cellCentre;     
  forAll(k.internalField(),celli)
    {
      cellCentre=mesh.C()[celli];
      z=cellCentre.z();
      k.internalField()[celli]=Foam::pow(ustar,2)/Foam::sqrt(Cmy);
      if(z-d+z0>0)
	epsilon.internalField()[celli]=Foam::pow(ustar,3)/(kappa*(z-d+z0));
      else
	epsilon.internalField()[celli]=Foam::pow(ustar,3)/(kappa*(z0));
      if ( (z-d)>0 && z<top )
	U[celli] =  U10*Foam::log((z-d+z0)/z0)/Foam::log((10.+z0)/z0);
      else if (z-d<0)
	U[celli] = U10*0.0;
      else
	U[celli] =  U10*Foam::log((top-d+z0)/z0)/Foam::log((10.+z0)/z0);
    }    

    U.write();
    k.write();
    epsilon.write();
    
    Info<< "End\n" << endl;
    
    return(0);
}


// ************************************************************************* //
 
