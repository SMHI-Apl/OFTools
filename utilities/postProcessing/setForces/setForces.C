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
  patchForces

  Description


  \*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OSspecific.H"
#include <fstream>
#include "wallFvPatch.H"
#include "argList.H"
#include "dimensionedScalar.H"
#include "faceSet.H"
#include "boundBox.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "RASModel.H"
//#include "incompressible/RASModel/RASModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  //argList::validArgs.append(".asc file");
  //argList::validOptions.insert("raster", "raster file");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readSetForcesDict.H"
#   include "createNu.H"
  
    dimensionedScalar rho("rho",dimensionSet(1,-3,0,0,0,0,0),1.25);
    dimensionedScalar mu(nu*rho);
  
    Info<<"Calculating forces for face sets"<<endl;
  
  
    const fvPatchList& patches=mesh.boundary();
  
    PtrList<faceSet> sets(setNames.size()); 
  
    forAll(setNames,setNamei)
    {
        faceSet fSet
        (
            mesh,
            setNames[setNamei],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        sets.set
        (
            setNamei, 
            new faceSet
            (
                mesh,
                setNames[setNamei],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    }
  
    List<vector> viscousForces(setNames.size());
    List<vector> pressureForces(setNames.size());
    List<scalar> areas(setNames.size());
  
    forAll(setNames,seti)
    {
        viscousForces[seti]=vector::zero;
        pressureForces[seti]=vector::zero;
        areas[seti]=0;
    }
  

    const surfaceVectorField::GeometricBoundaryField& Sfb =
          mesh.Sf().boundaryField();

    tmp<volSymmTensorField> tdevRhoReff = rho*RASModel->devReff();
    const volSymmTensorField::GeometricBoundaryField& devRhoReffb = tdevRhoReff().boundaryField();

    forAll(sets,seti)
    {
       Info<<"Processing set: "<<setNames[seti]<<endl;
       label foundInPatch=0;
       forAll(patches,patchi)
       {
	 Info<<"Searching for cell faces in patch: "<<patches[patchi].name()<<endl;
           
           const fvPatch& curPatch=patches[patchi];
           const polyPatch& pp=mesh.boundaryMesh()[patchi];
           const pointField &centres=mesh.C().boundaryField()[patchi];
           
           scalarField mask(centres.size(),0);         
           forAll(curPatch.Sf(),i)
           {   

               label globInd=i+pp.start();
               for(faceSet::iterator facei=sets[seti].begin();facei!=sets[seti].end();facei++)
               {
                   label faceInd=facei.key();
                   if(faceInd==globInd)
                   {
                       mask[i]=1;
		       foundInPatch+=1;
                   }
               }
           }
  
	   
           scalar area;
           vector pressureForce;
           vector viscousForce;
  
           pressureForce = sum
           (
               mask*p.boundaryField()[patchi]*rho.value() 
             * Sfb[patchi]
           );
      
           viscousForce = sum(mask*devRhoReffb[patchi] & Sfb[patchi]);
           
           area=sum
           (
              mask*mesh.magSf().boundaryField()[patchi]
           );
           
           pressureForces[seti]=pressureForces[seti]+pressureForce;
           
           viscousForces[seti]=viscousForces[seti]+viscousForce;
           
           areas[seti]=areas[seti]+area;
       }
    }
  
  
    Info<<"Set\tarea\tFv_x[N]\tFv_y[N]\tFv_z[N]\tFp_x[N]\tFp_y[N]\tFp_z[N]"<<endl;
    forAll(setNames,seti)
    {
        Info << setNames[seti]<<"\t"
             << areas[seti]<<"\t"
             << viscousForces[seti][0]<<"\t"
             << viscousForces[seti][1]<<"\t"
             << viscousForces[seti][2]<<"\t"
             << pressureForces[seti][0]<<"\t"
             << pressureForces[seti][1]<<"\t"
             << pressureForces[seti][2]<<"\t"
             << endl;
    }      
  
    Info<<"\nCalculating forces for face sets finished succesfully!\n"<<endl;
    return(0);
}

// ************************************************************************* //

