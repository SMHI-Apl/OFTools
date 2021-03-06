/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::incompressible::RASModels::kEpsilon_canopy_veh

Group
    grpIcoRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible flows.

    Reference:
    \verbatim
        Launder, B.E., and Spalding, D.B.,
        "Mathematical Models of Turbulence",
        Academic Press, 1972.

        Launder, B.E., and Spalding, D.B.,
        "Computational Methods for Turbulent Flows",
        Comp. Meth in Appl Mech & Eng'g, Vol 3, 1974, pp 269-289.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        kEpsilon_canopy_vehCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    kEpsilon_canopy_veh.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilon_canopy_veh_H
#define kEpsilon_canopy_veh_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilon_canopy_veh Declaration
\*---------------------------------------------------------------------------*/

class kEpsilon_canopy_veh
:
    public RASModel
{

protected:

    // Protected data

        // Coefficients for Traffic Produced Turbulence
            dimensionedScalar Cfcar_;
            dimensionedScalar L_;
            dimensionedScalar Cecar_;

        // Canopy coefficients
            dimensionedScalar betaP_;
            dimensionedScalar betaD_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;


        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaEps_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;

        // Fields representing vegetation canopy
    
            volScalarField canopy_;

        //Fields representing vehicle canopy

            volVectorField Ucar_;
            volScalarField VAD_;

public:

    //- Runtime type information
    TypeName("kEpsilon_canopy_veh");

    // Constructors

        //- Construct from components
        kEpsilon_canopy_veh
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kEpsilon_canopy_veh()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the canopy field
        volScalarField& canopy()
        {
            return canopy_;
        }
    
        //- Return the car velocity field
        virtual tmp<volVectorField> Ucar() const
        {
            return Ucar_;
        }

        //- Return the Vehicle Area Density field
        virtual tmp<volScalarField> VAD() const
        {
            return VAD_;
        }


        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }


        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
