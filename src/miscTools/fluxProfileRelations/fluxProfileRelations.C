/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fluxProfileRelations.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxProfileRelations, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxProfileRelations::fluxProfileRelations
(
    word relTypeName
)
:
relType(SMOOTH)

{
    if (relTypeName == "smooth" ||
        relTypeName == "Smooth" ||
        relTypeName == "smoothWall" ||
        relTypeName == "SmoothWall" ||
        relTypeName == "smoothwall")
    {
        relType = SMOOTH;
    }
    else if (relTypeName == "moninobukhov" ||
             relTypeName == "Moninobukhov" ||
             relTypeName == "MoninObukhov" ||
             relTypeName == "MO" ||
             relTypeName == "mo")
    {
        relType = MONIN_OBUKHOV;
    }

  //Info << "relType = " << relType << endl;
}




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxProfileRelations::~fluxProfileRelations()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<scalar> Foam::fluxProfileRelations::update
(
    const scalar zref_,
    const scalar Uref_,
    const scalar kappa_,
    const scalar B_,
    const scalar nu_,
    const scalar beta_m_,
    const scalar beta_h_,
    const scalar gamma_m_,
    const scalar gamma_h_,
    const scalar alpha_h_,
    const scalar z0_,
    const scalar Tref_,
    const scalar Tsurf_,
    const scalar Qsurf_,
    const word surfTempOrFluxName
)
{
    // Update the class variables using those passed into this function.
    zref = zref_;
    Uref = Uref_;
    Tref = Tref_;
    Tsurf = Tsurf_;
    Qsurf = Qsurf_;
    kappa = kappa_;
    B = B_;
    nu = nu_;
    beta_m = beta_m_;
    beta_h = beta_h_;
    gamma_m = gamma_m_;
    gamma_h = gamma_h_;
    alpha_h = alpha_h_;
    z0 = z0_;


    // Set the surface temperature or flux enumerator.
    surfaceTempOrFlux surfTempOrFlux = TEMPERATURE;
 
    if (surfTempOrFluxName == "temperature" ||
        surfTempOrFluxName == "Temperature" ||
        surfTempOrFluxName == "temp" ||
        surfTempOrFluxName == "Temp")
    {
        surfTempOrFlux = TEMPERATURE;
    }
    else if (surfTempOrFluxName == "flux" ||
             surfTempOrFluxName == "Flux" ||
             surfTempOrFluxName == "tempFlux" ||
             surfTempOrFluxName == "TempFlux" ||
             surfTempOrFluxName == "temperatureFlux" ||
             surfTempOrFluxName == "TemperatureFlux")
    {
        surfTempOrFlux = FLUX;
    }


    // Initialize the list that will return friction velocity, surface temperature flux, and surface temperature.
    List<scalar> fluxes(3, 0.0);

 
    // Select and call the specific update type.
    switch (relType)
    {
        // Smooth wall case.
        case SMOOTH:
            fluxes = updateSmooth();
            break;

        // Rough wall Monin-Obukhov case.
        case MONIN_OBUKHOV:
            fluxes = updateMoninObukhov();
            break;

        // Default case.
        default:
            fluxes = updateMoninObukhov();
            break;
    }

    return fluxes;
}



List<scalar> Foam::fluxProfileRelations::updateSmooth()
{
    // We use a Newton-Raphson solver to solve the log-law equation that is nonlinear in 
    // terms of friction velocity.  This solver method was developed by Prakash Mohan (NREL).

    // Set the solver tolerance and maximum number of iterations.  If either the residual drops
    // below the tolerance, or the iterations exceeds the maximum, the solver exits.
    const scalar tol = 1.0E-5;
    const int iterMax = 20;
    scalar residual = 1.0;
    int iter = 0;
  
    // Set the initial guess for friction velocity.
    scalar utau = 0.1;
  //Info << "nu, B, kappa = " << nu << ", " << B << ", " << kappa << endl;

    // The solver loop.
    while ((Foam::mag(residual) > tol) && (iter < iterMax))
    {   
      //Info << "iter = " << iter << tab << "utau = " << utau << tab << "residual = " << residual << endl;

        // This is the derivative of velocity magnitude with respect to friction velocity.
        scalar fPrime = ((1.0/kappa)*(1.0 + Foam::log(utau*zref/nu)) + B);

        // The residual used to update the friction velocity.
        residual = (utau*((1.0/kappa)*Foam::log(utau*zref/nu) + B) - Uref) / fPrime;

        // Update the friction velocity.
        utau -= residual;
    }
 
    // Return the friction velocity, surface temperature flux, and surface temperature.   
    List<scalar> fluxes(3, 0.0);
    fluxes[0] = utau;
    fluxes[1] = Qsurf;
    fluxes[2] = Tsurf;

    return fluxes;
}

List<scalar> Foam::fluxProfileRelations::updateMoninObukhov()
{
    List<scalar> fluxes(3, 0.0);
    return fluxes;
}

// ************************************************************************* //
