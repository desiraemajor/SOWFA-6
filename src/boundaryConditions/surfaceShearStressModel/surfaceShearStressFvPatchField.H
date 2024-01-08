/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Class
    Foam::incompressible::LESModels::surfaceShearStressFvPatchField

Description
    Wall boundary condition for total shear stress at a solid boundary.  

SourceFiles
    surfaceShearStressFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceShearStressFvPatchField_H
#define surfaceShearStressFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedFvPatch.H"
//#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class surfaceShearStressFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class surfaceShearStressFvPatchField
:
    public fixedValueFvPatchSymmTensorField//,
//    public mappedPatchBaseVectorFieldField
{

private:

    // Private data
        enum class averagingType
        {
            NONE,
            PLANAR,
            TEMPORAL
        };

        enum class fluctuationModel
        {
            NONE,
            SCHUMANN,
            MOENG
        };

        //- Various flux profile relation constants.
        scalar kappa;

        scalar beta_m;
        scalar beta_h;
        scalar gamma_m;
        scalar gamma_h;
        scalar alpha_h;

        scalar B;

        //- Roughness height (which can be space varying).
        scalarField z0;

        //- Type of averaging used on velocity/temperature before passing to flux 
        //  profile relations.
        //  Options are:
        //    --none
        //    --planar
        //    --temporal
        word avgTypeName;
        averagingType avgType;

        //- Model to be applied to make stress field fluctuating.
        word fluctModelName;
        fluctuationModel fluctModel;

        //- Molecular viscosity.
        scalar nu;

        //- Interpolation scheme used for 'nearestCell' sampling.
        word interpolationScheme;

        vectorField USampled;



    // Private member functions
       
/*
protected:

    // Protected Member Functions

        const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<vector, volMesh>& iF
        );
*/

public:

    //- Runtime type information
    TypeName("surfaceShearStressModel");


    // Constructors

        //- Construct from patch and internal field
        surfaceShearStressFvPatchField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceShearStressFvPatchField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  stressABLWallFunctionFvPatchScalarField
        //  onto a new patch
        surfaceShearStressFvPatchField
        (
            const surfaceShearStressFvPatchField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceShearStressFvPatchField
        (
            const surfaceShearStressFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new surfaceShearStressFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceShearStressFvPatchField
        (
            const surfaceShearStressFvPatchField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new surfaceShearStressFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access functions
            const scalarField& z0_() const
            {
                return z0;
            }
            scalarField& z0_()
            {
                return z0;
            }

        // Evaluation functions

            void sampleVelocityField(vectorField& USampled);

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes); 
          //virtual void evaluate();
          //virtual void updateCoeffs();

        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
