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
    transformTopography

Description
    Transforms the mesh points in the polyMesh directory according to the
    given topography:

    The height to follow ground strictly is given as an option. Above this height the mesh is graded in order to get a "flat" upper boundary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "transformField.H"
#include "IStringStream.H"
#include "Raster.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append(".asc file");
  //  argList::validOptions.insert("height", "scalar");
    
#   include "setRootCase.H"
#   include "createTime.H"
#   include "readLanduseDict.H"
    fileName topoName(args.args()[1]);

    Raster topo;

    if(!topo.read(topoName.c_str()))
    FatalErrorIn(args.executable())<< "Cannot read file "<< topoName << exit(FatalError);

    topo.xll=topo.xll-double(subtractedX);
    topo.xur=topo.xur-double(subtractedX);
    topo.yll=topo.yll-double(subtractedY);
    topo.yur=topo.yur-double(subtractedY);


    Info<<"-----------Topography specification-------------"<<endl;
    Info<<"Extent: "<<topo.xll<<"< X <"<<topo.xur<<" and "
	<<topo.yll<<"< Y <"<<topo.yur<<endl;
    
    Info<<"Dimensions: "<<"cellsize="<<topo.cellsize
	<<" nrows= "<<topo.nrows<<" ncols= "<<topo.ncols<<endl; 


    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    scalar roof=-9999;
    scalar offset=9e99;
    scalar offsetMax=-9e99;
    forAll(points,pInd)
      {
	vector& p=points[pInd];
	if(p.z() > roof)
	  roof=p.z();
	
	scalar ground = topo.interpValue(double(p.x()),double(p.y()));
	if(ground < offset)
	  offset=ground;
	if(ground > offsetMax)
	  offsetMax=ground;
      }

    Info<< "The maximum height of the mesh is "<< roof << endl;
    Info<< "The min height of the topography in the domain is "<< offset << endl;
    Info<< "The max height of the topography in the domain is "<< offsetMax << endl;
    
    roof = roof + offset;

    scalar dampedMove0=-9999;

    forAll(points,pInd)
      {
	vector& p=points[pInd];
	scalar ground       = topo.interpValue(double(p.x()),double(p.y()));

	// move all points to lowest level of topo
	p.z() = p.z() + offset;
	
	// points will be moved further by maximum (ground-offset)
	scalar maxMove = ground - offset;
        scalar dampedMove   = maxMove * (roof-p.z())/(roof-offset);
        
	p.z() = p.z() + dampedMove;

	if(dampedMove > dampedMove0)
	  dampedMove0=dampedMove;


      }

    IOstream::defaultPrecision(10);

    Info << "The maximum movement of mesh is "<< dampedMove0 << endl;
    Info << "Writing modified points into directory " << points.path() << nl << endl;
    points.write();

    return(0);
}


// ************************************************************************* //

