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

    topo.xll = topo.xll - double(subtractedX);
    topo.xur = topo.xur - double(subtractedX);
    topo.yll = topo.yll - double(subtractedY);
    topo.yur = topo.yur - double(subtractedY);


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

   // scalar height=50; not in use

  /* if (args.options().found("height"))
      {
	IStringStream instr(args.options()["height"]);
	instr >> height;
      }
    Info<< "Points are set to follow ground up to " << height << " m"<<endl;
 */
    scalar roof=-9999;
    scalar offset=9e99;
    forAll(points,pInd)
      {
	vector& p=points[pInd];
	if(p.z() > roof)
	  roof=p.z();
	
	scalar ground = topo.interpValue(double(p.x()),double(p.y()));
	if(ground < offset)
	  offset=ground;
      }

    Info<< "The maximum height of the mesh is "<< roof << endl;
    Info<< "The min height of the topography in the domain is "<< offset << endl;
    
    roof = roof + offset;
    forAll(points,pInd)
      {
	vector& p = points[pInd];
	scalar ground       = topo.interpValue(double(p.x()), double(p.y()));
        scalar undampedMove = ground - offset;
        scalar dampedMove   = undampedMove * (roof - p.z()) / (roof - ground);
        
	p.z() = p.z() + offset;
	p.z() = p.z() + dampedMove;
	// p.z()=(p.z()-height)*(roof-ground-height)/(roof-height)+ground+height;
      }

    //    if (args.options().found("scale"))
    //{
    //    vector scaleVector(IStringStream(args.options()["scale"])());
    //
    //  Info<< "Scaling points by " << scaleVector << endl;
    //
    //  points.replace(vector::X, scaleVector.x()*points.component(vector::X));
    //  points.replace(vector::Y, scaleVector.y()*points.component(vector::Y));
    //  points.replace(vector::Z, scaleVector.z()*points.component(vector::Z));
    //}

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing modified points into directory " << points.path() << nl << endl;
    points.write();

    return(0);
}


// ************************************************************************* //

