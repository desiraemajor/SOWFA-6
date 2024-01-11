/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fluxProfileRelations

Description
    A parent class for flux profile relations, such as Monin-Obukhov 
    similarity.

SourceFiles
    fluxProfileRelations.C

\*---------------------------------------------------------------------------*/

#ifndef fluxProfileRelations_H
#define fluxProfileRelations_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class meshPlanes Declaration
\*---------------------------------------------------------------------------*/

class fluxProfileRelations
{
    // Private data
    // - An enumeration holding the type of flux relation to use.
    enum relationType
    {
        SMOOTH,
        MONIN_OBUKHOV,
        DONELAN
    };

    // - An enumeration containing the user's selection of supplying
    //   surface temperature or surface temperature flux.
    enum surfaceTempOrFlux
    {
        FLUX,
        TEMPERATURE
    };

    // - Type of flux profile relation
    relationType relType;

    // - Flux relation constants.
    // - von Karman constant
    scalar kappa;

    // - constant in smooth wall log law
    scalar B;

    // - molecular viscosity
    scalar nu;

    // - constants in Monin-Obukhov relations
    scalar beta_m;
    scalar beta_h;
    scalar gamma_m;
    scalar gamma_h;
    scalar alpha_h;
    scalar z0;

    // - Flux relation velocities, temperatures, fluxes, sampling height
    scalar zref;
    scalar Uref;
    scalar Tref;
    scalar Tsurf;
    scalar Qsurf;
    word surfTempOrFluxName;



    // Private Member Functions



public:

    //- Declare name of the class and its debug switch
    ClassName("fluxProfileRelations");




    // Static member data




    // Constructors

        //- Construct from components
        fluxProfileRelations
        (
            word relTypeName
        );




    // Destructor
    virtual ~fluxProfileRelations();




    // Member functions

    // - This is the overall update function that calls specific update functions
    //   depending on the type of flux profile relationship the user has asked
    //   for, which is stored in "relType".  It returns a triple list of scalars
    //   containing friction velocity, surface temperature flux, and surface
    //   temperature.  Not all relations compute surface temperature flux and
    //   surface temperature, so in that case, the zero flux values are returned.
    List<scalar> update
    (
        const scalar zref_,
        const scalar Uref_,
        const scalar kappa_ = 0.4,
        const scalar B_ = 5.0,
        const scalar nu_ = 1.0E-5,
        const scalar beta_m_ = 16.0,
        const scalar beta_h_ = 16.0,
        const scalar gamma_m_ = 5.0,
        const scalar gamma_h_ = 5.0,
        const scalar alpha_h_ = 0.9,
        const scalar z0_ = 0.1,
        const scalar Tref_ = 300.0, 
        const scalar Tsurf_ = 300.0,
        const scalar Qsurf_ = 0.0,
        const word surfTempOrFluxName = "temperature"
    );

    // - This is the update function for smooth walls that uses the smooth-
    //   wall log profile to give friction velocity given a velocity 
    //   magnitude at a reference height on the log profile.
    List<scalar> updateSmooth();

    // - This is the update function for rough wall atmospheric flows that
    //   uses the similarity laws of Monin and Obukhov to return friction
    //   velocity and surface temperature flux or surface temperature given the 
    //   velocity magnitude and potential temperature at reference height or
    //   surface heat flux. 
    List<scalar> updateMoninObukhov();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
