/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    liftDrag

Description
    Calculates the pressure loads on the Kista Torns. Patch names,
    centres of moments etc. are all hard coded.

    Note the possibility to define bounding boxes in order to only include
    patch areas inside these boxes in the calculations, e.g. if only the 10
    upper meters of one tower is of interest.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dimensionedTypes.H"
#include "incompressible/transportModel/transportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/LESmodel/LESmodel.H"

#include "wallFvPatch.H"
#include "liftDrag.H"

scalar xyComponent(const vector v)
{
    return Foam::sqrt(Foam::pow(v.x(),2.0)+Foam::pow(v.y(),2.0));
}

scalar xyAngle(const vector v)
{
    scalar toDeg = 180/M_PI;
    scalar angle = toDeg*Foam::atan(v.y()/(v.x()+SMALL));
    
    if ( (v.y() > 0) && (v.x() < 0) ) angle=360-angle; 
    if ( (v.y() < 0) && (v.x() < 0) ) angle+=180; 
    if ( (v.y() < 0) && (v.x() > 0) ) angle+=90; 
    return angle;
}


pointField getPatchPoints(const fvMesh &mesh, const word &patchName)
{

    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const labelList &patchPointLabels(mesh.boundaryMesh()[patchID].meshPoints());
    pointField tmpPoints(patchPointLabels.size());
    forAll (patchPointLabels, i )
    {
        tmpPoints[i] = mesh.points()[patchPointLabels[i]];
    }
    return tmpPoints;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"


    /*
     * Configuration starts here
     */
    wordList patchNames(2);
    patchNames[0]="High_tower_k1";
    patchNames[1]="Lower_tower_k2";

    // Points around which to calculate moments
    // one for each patch (tower)
    List<vector> cPoints(2);
    cPoints[0]= vector(13.2475,-23.525,18);
    cPoints[1]= vector(-19.925,-53.930,18);

    scalar deltaZ=1;
    scalar min=-1000;
    scalar max=1000;
    scalar densityCorrFactor=1.25;
    scalar velocityCorrFactor=0.4489;
    scalar offsetHeight=18;

    /*
     * End configuration
     */


#   include "createNu.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject pHeader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check p exists
        if (pHeader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading p\n" << endl;
            volScalarField p(pHeader, mesh);

            forAll(patchNames, patchNameI)
            {
                word patchName=patchNames[patchNameI];
                vector cPoint=cPoints[patchNameI];

                boundBox patchBox(getPatchPoints(mesh,patchName));

                scalar Zmin = patchBox.min().z();
                scalar Zmax = patchBox.max().z();
                scalar ZcurMin=Zmin;
                scalar ZcurMax=ZcurMin+deltaZ;

                Info<< "\nWall patch " << " named " << patchName << " : " << nl
                    << "Vertical extent [m]: " << Zmin << " to " << Zmax << nl
                    << "Vertical integration span width [m]: " << deltaZ << endl;
                while (ZcurMin <= Zmax)
                {
                    boundBox curBox(vector(min,min,ZcurMin),vector(max,max,ZcurMax));
                    vector pl = pressureLoad(p,patchName,curBox)*velocityCorrFactor*densityCorrFactor;
                    Info << "From " << ZcurMin-offsetHeight << " to " << ZcurMax-offsetHeight << ": "<< pl << endl;

                    // Climbing one step:
                    ZcurMin=ZcurMax;
                    ZcurMax+=deltaZ;
                }
            }
        }
        else
        {
            Info<< "    No p" << endl;
        }

        Info<< endl;
    }

    return(0);
}


// ************************************************************************* //
