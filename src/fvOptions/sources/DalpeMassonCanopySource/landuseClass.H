/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
    
    Copyright of this code belong to the author.    

Class
    Foam::landuseClass

Author
    David Segersson (david.segersson@smhi.se)
    Swedish Meteorological and Hydrological Institute

Description
    A landuseClass

SourceFiles
    landuseClass.H
    landuseClass.C

\*---------------------------------------------------------------------------*/

#ifndef landuseClass_H
#define landuseClass_H

#include "label.H"
#include "scalar.H"
#include "dimensionedScalar.H"
#include "dimensionSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class landuseClass Declaration
\*---------------------------------------------------------------------------*/

class landuseClass {

  private:
    label code_;
    word name_;
    scalar Cd_;
    scalar LAI_;
    scalar z0_;
    scalar height_;
    scalar LADmax_;
    scalarList LADProfile_;
  
  public:
  
    // Constructors
    landuseClass();
    landuseClass(const dictionary &dict, word name);

    // Destructor
    ~landuseClass();

    // member functions
    const scalarList& LADProfile() {
      return LADProfile_;
    }

    void LADmaxFromLAI();
    scalar integrateLAD();

    // access
    const word& name() {
      return name_;
    }

    const label& code() {
      return code_;
    }

    const scalar& Cd() {
      return Cd_;
    }

    const scalar& LAI() {
      return LAI_;
    }

    const scalar& z0() {
      return z0_;
    }

    const scalar& height() {
      return height_;
    }

    const scalar& LADmax() {
      return LADmax_;
    }  

    
};
}

#endif /*landuseClass_H*/
