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

\*---------------------------------------------------------------------------*/

#include "surfaceShearStressFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "geometricTransformations.H"
#include "fluxProfileRelations.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceShearStressFvPatchField::
surfaceShearStressFvPatchField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(p, iF),
//  mappedPatchFieldBase<vector>(this->mapper(p, iF), *this),
    kappa(0.41),
    beta_m(16.0),
    beta_h(16.0),
    gamma_m(5.0),
    gamma_h(5.0),
    alpha_h(0.9),
    B(5.0),
    z0(p.size(), 0.1),
    fluxProfileRelTypeName("MoninObukhov"),
    avgTypeName("none"),
    avgType(averagingType::NONE),
    fluctModelName("none"),
    fluctModel(fluctuationModel::NONE),
    nu(8.0E-6),
    interpolationScheme("cell"),
    patchName("lower")
{
    Info << "In constructor of surfaceShearStressFvPatchField (1)..." << endl;
    USampled = vectorField(p.size(), Zero);
}


Foam::surfaceShearStressFvPatchField::
surfaceShearStressFvPatchField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(p, iF, dict),
//  mappedPatchFieldBase<vector>(this->mapper(p, iF), *this, ptf),
    kappa(dict.lookupOrDefault<scalar>("kappa",0.41)),
    beta_m(dict.lookupOrDefault<scalar>("beta_m",16.0)),
    beta_h(dict.lookupOrDefault<scalar>("beta_h",16.0)),
    gamma_m(dict.lookupOrDefault<scalar>("gamma_m",5.0)),
    gamma_h(dict.lookupOrDefault<scalar>("gamma_h",5.0)),
    alpha_h(dict.lookupOrDefault<scalar>("alpha_h",0.9)),
    B(dict.lookupOrDefault<scalar>("B",5.0)),
  //z0(dict.lookupOrDefault<scalarField>("z0",scalarField(p.size(),0.1))),
  //z0("z0", dict, p.size()),
    fluxProfileRelTypeName(dict.lookupOrDefault<word>("fluxProfileRelationType","MoninObukhov")),
    avgTypeName(dict.lookupOrDefault<word>("averageType","none")),
    avgType(averagingType::NONE),
    fluctModelName(dict.lookupOrDefault<word>("fluctuationModel","none")),
    fluctModel(fluctuationModel::NONE),
    nu(8.0E-6),
    interpolationScheme(dict.lookupOrDefault<word>("interpolationScheme","cell"))
{
    if (dict.found("z0"))
    {
        z0 = scalarField("z0", dict, p.size());
    }
    else
    {
        z0 = scalarField(p.size(), 0.1);
    }

    patchName = this->patch().name();

    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
    Info << "In constructor of surfaceShearStressFvPatchField (2)..." << endl;
  //Info << "z0 = " << z0 << endl;
  //Info << "z0 = " << z0_() << endl;
  //sampleVelocityField(USampled);
    USampled = vectorField(p.size(), Zero);
  //Info << "USampled.size() = " << USampled.size() << endl;

    openOutputFile();
}


Foam::surfaceShearStressFvPatchField::
surfaceShearStressFvPatchField
(
    const surfaceShearStressFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(ptf, p, iF, mapper),
//  mappedPatchFieldBase<vector>(this->mapper(p, iF), *this, ptf),
    kappa(ptf.kappa),
    beta_m(ptf.beta_m),
    beta_h(ptf.beta_h),
    gamma_m(ptf.gamma_m),
    gamma_h(ptf.gamma_h),
    alpha_h(ptf.alpha_h),
    B(ptf.B),
    z0(ptf.z0, mapper),
    fluxProfileRelTypeName(ptf.fluxProfileRelTypeName),
    avgTypeName(ptf.avgTypeName),
    avgType(ptf.avgType),
    fluctModelName(ptf.fluctModelName),
    fluctModel(ptf.fluctModel),
    nu(ptf.nu),
    interpolationScheme(ptf.interpolationScheme),
    USampled(ptf.USampled),
    patchName(ptf.patchName)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
    Info << "In constructor of surfaceShearStressFvPatchField (3)..." << endl;
}


Foam::surfaceShearStressFvPatchField::
surfaceShearStressFvPatchField
(
    const surfaceShearStressFvPatchField& ptf
)
:
    fixedValueFvPatchSymmTensorField(ptf),
//  mappedPatchFieldBase<vector>(ptf),
    kappa(ptf.kappa),
    beta_m(ptf.beta_m),
    beta_h(ptf.beta_h),
    gamma_m(ptf.gamma_m),
    gamma_h(ptf.gamma_h),
    alpha_h(ptf.alpha_h),
    B(ptf.B),
    z0(ptf.z0),
    fluxProfileRelTypeName(ptf.fluxProfileRelTypeName),
    avgTypeName(ptf.avgTypeName),
    avgType(ptf.avgType),
    fluctModelName(ptf.fluctModelName),
    fluctModel(ptf.fluctModel),
    nu(ptf.nu),
    interpolationScheme(ptf.interpolationScheme),
    USampled(ptf.USampled),
    patchName(ptf.patchName)
{
    Info << "In constructor of surfaceShearStressFvPatchField (4)..." << endl;
}


Foam::surfaceShearStressFvPatchField::
surfaceShearStressFvPatchField
(
    const surfaceShearStressFvPatchField& ptf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(ptf, iF),
//  mappedPatchFieldBase<vector>(this->mapper(this->patch(), iF), *this, ptf),
    kappa(ptf.kappa),
    beta_m(ptf.beta_m),
    beta_h(ptf.beta_h),
    gamma_m(ptf.gamma_m),
    gamma_h(ptf.gamma_h),
    alpha_h(ptf.alpha_h),
    B(ptf.B),
    z0(ptf.z0),
    fluxProfileRelTypeName(ptf.fluxProfileRelTypeName),
    avgTypeName(ptf.avgTypeName),
    avgType(ptf.avgType),
    fluctModelName(ptf.fluctModelName),
    fluctModel(ptf.fluctModel),
    nu(ptf.nu),
    interpolationScheme(ptf.interpolationScheme),
    USampled(ptf.USampled),
    patchName(ptf.patchName)
{
    Info << "In constructor of surfaceShearStressFvPatchField (5)..." << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
void Foam::surfaceShearStressFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->mappedField());

    if (debug)
    {
        Info<< "mapped on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchSymmTensorField::updateCoeffs();
}
*/



//void Foam::surfaceShearStressFvPatchField::updateCoeffs()
//void Foam::surfaceShearStressFvPatchField::evaluate()
void Foam::surfaceShearStressFvPatchField::evaluate(const Pstream::commsTypes)
{
/*
    if (updated())
    {
        return;
    }
*/
    sampleVelocityField(USampled);

    const vectorField normal = -patch().nf();
    const scalarField area = patch().magSf();
    const scalar areaSum = gSum(area);

    symmTensorField& Rw = *this;

  //Info << "Rw.size() = " << Rw.size() << endl;
  //Info << "USampled.size() = " << USampled.size() << endl;

    scalar utauMean = 0.0;

    const dictionary& transportProperties = db().lookupObject<dictionary>("transportProperties");
    dimensionedScalar nuDim = transportProperties.lookup("nu");
    nu = nuDim.value();

    
    forAll(Rw, faceI)
    {
        scalar u = USampled[faceI].x();
        scalar v = USampled[faceI].y();
        scalar U = max(Foam::sqrt(Foam::sqr(u) + Foam::sqr(v)), 1.0E-5);

        // Declare an instance of the flux-profile-relations class.
        fluxProfileRelations fluxProfileRel(fluxProfileRelTypeName);

        // Get the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            surfaceShearStressFvPatchField::patch().patch()
        );
        vector offset = mpp.offset();
      //Info << "offset = " << offset << endl;

      //Info << "samplePoints() = " << mpp.samplePoints() << endl;

        scalar zRef = normal[faceI] & offset;
      //Info << "zRef = " << zRef << endl;
        List<scalar> fluxes = fluxProfileRel.update(zRef, U, kappa, B, nu);
        scalar utau = fluxes[0];
        utauMean += utau * area[faceI];

        Rw[faceI] = Zero;
        Rw[faceI].xz() = -normal[faceI].z() * Foam::sqr(utau) * (u / U);
        Rw[faceI].yz() = -normal[faceI].z() * Foam::sqr(utau) * (v / U);
    }

    reduce(utauMean,sumOp<scalar>());
    utauMean = (areaSum > 0.0) ? utauMean/areaSum: 0.0;

    Info << "<u_tau> = " << utauMean << " m/s" << endl;

    writeOutputFile(utauMean, 300.0, 0.0);
  //Info << USampled[0] << tab << endl;

  //fixedValueFvPatchSymmTensorField::updateCoeffs();
}



//void Foam::surfaceShearStressFvPatchField::updateCoeffs()
void Foam::surfaceShearStressFvPatchField::sampleVelocityField(vectorField& USampled)
{
  //if (updated())
  //{
  //    return;
  //}

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;


    // Get the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        surfaceShearStressFvPatchField::patch().patch()
    );
    const fvMesh& thisMesh = surfaceShearStressFvPatchField::patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
  //const word& fieldName = internalField().name();
  //const volVectorField& UField =
  //    nbrMesh.lookupObject<volVectorField>(fieldName);
    const volVectorField& UField =
        nbrMesh.lookupObject<volVectorField>("U");

  //word phiName_ = "phi";
  //surfaceScalarField& phiField =
  //    nbrMesh.lookupObjectRef<surfaceScalarField>(phiName_);

  //vectorField USampled;
  //scalarField newPhiValues;
  //symmTensorField newRwValues;
  //symmTensorField& Rw = *this;

    vectorField samples(mpp.samplePoints());

  //newRwValues.setSize(samples.size(), pTraits<symmTensor>::zero);

    switch (mpp.mode())
    {
        case mappedPolyPatch::NEARESTCELL:
        {
            vectorField allUValues(nbrMesh.nFaces(), Zero);

            const mapDistribute& distMap = mpp.map();

          //Info << "distMap = " << distMap << endl;
          //Info << "nbrMesh.nFaces() = " << nbrMesh.nFaces() << endl;
          //Info << "allUValues = " << allUValues << endl;
          //Info << "USampled = " << newUValues << endl;
          //Info << "newRwValues = " << newRwValues << endl;

            if (interpolationScheme != interpolationCell<vector>::typeName)
            {
              //Info << "In if (1)..." << endl;
                // Send back sample points to the processor that holds the cell
              //vectorField samples(mpp.samplePoints());
                distMap.reverseDistribute
                (
                    (
                        mpp.sameRegion()
                      ? thisMesh.nCells()
                      : nbrMesh.nCells()
                    ),
                    point::max,
                    samples
                );

                autoPtr<interpolation<vector>> interpolator
                (
                    interpolation<vector>::New
                    (
                        interpolationScheme,
                        UField
                      //sampleField()
                    )
                );
                const interpolation<vector>& interp = interpolator();
                
                USampled.setSize(samples.size(), pTraits<vector>::max);
              //newRwValues.setSize(samples.size(), pTraits<symmTensor::zero);
                forAll(samples, celli)
                {
                    if (samples[celli] != point::max)
                    {
                        allUValues[celli] = interp.interpolate
                        (
                            samples[celli],
                            celli
                        );
                    }
                }

            }
            else
            {
              //Info << "In if (2)..." << endl;


              //typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
 
              //const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
 
                vectorField samples(mpp.samplePoints());
              //newRwValues.setSize(samples.size(), pTraits<symmTensor>::max);
              //newRwValues.setSize(samples.size(), pTraits<symmTensor>::zero);
                if (mpp.sameRegion())
                {
                  //Info << "In mpp.sameRegion()..." << endl;
                    /*
                    if (fieldName_ == patchField_.internalField().name())
                    {
                        // Optimisation: bypass field lookup
                        return
                            dynamic_cast<const fieldType&>
                            (
                                patchField_.internalField()
                            );
                    }
                    */
                  //else
                  //{
                      //const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
                      //return thisMesh.template lookupObject<fieldType>(fieldName_);
                    allUValues = thisMesh.lookupObject<volVectorField>("U");
                  //Info << "allUValues = " << allUValues << endl;
                  //}
                }
                else
                {
                  //Info << "Not in mpp.sameRegion()..." << endl;
                  //return nbrMesh.template lookupObject<fieldType>(fieldName_);
                    allUValues = nbrMesh.lookupObject<volVectorField>("U");
                }


              //allUValues = sampleField();
            }

            distMap.distribute(allUValues);
            USampled.transfer(allUValues);

          //Info << "USampled = " << newUValues << endl;
          //Info << "newRwValues = " << newRwValues << endl;

            break;
        }
        case mappedPolyPatch::NEARESTFACE:
        {
            vectorField allUValues(nbrMesh.nFaces(), Zero);
          //scalarField allPhiValues(nbrMesh.nFaces(), 0.0);
          //symmTensorField allRwValues(nbrMesh.nFaces(), Zero);

            forAll(UField.boundaryField(), patchi)
            {
                const fvPatchVectorField& Upf = UField.boundaryField()[patchi];
              //const scalarField& phipf = phiField.boundaryField()[patchi];

                label faceStart = Upf.patch().start();

                forAll(Upf, facei)
                {
                    allUValues[faceStart + facei] = Upf[facei];
                  //allPhiValues[faceStart + facei] = phipf[facei];
                }
            }

            mpp.distribute(allUValues);
            USampled.transfer(allUValues);

          //mpp.distribute(allPhiValues);
          //newPhiValues.transfer(allPhiValues);

            break;
        }
        case mappedPolyPatch::NEARESTPATCHFACE:
        case mappedPolyPatch::NEARESTPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mpp.samplePatch());

            USampled = UField.boundaryField()[nbrPatchID];
            mpp.distribute(USampled);

          //newPhiValues = phiField.boundaryField()[nbrPatchID];
          //mpp.distribute(newPhiValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "patch can only be used in NEARESTPATCHFACE, "
                << "NEARESTPATCHFACEAMI or NEARESTFACE mode" << nl
                << abort(FatalError);
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::surfaceShearStressFvPatchField::openOutputFile()
{
    const fvMesh& thisMesh = surfaceShearStressFvPatchField::patch().boundaryMesh().mesh();

    fileName outputPath(fileName::null);

    // Specify output path
    if (Pstream::parRun())
    {
        outputPath = thisMesh.time().path()/".."/"postProcessing"/("surfaceFluxes." & patchName)/thisMesh.time().timeName();
    }
    else
    {
        outputPath = thisMesh.time().path()/"postProcessing"/("surfaceFluxes." & patchName)/thisMesh.time().timeName();
    }

    if (Pstream::master)
    {
        // Create directory if it does not exist
        mkDir(outputPath);

        word fileName;

        fileName = "surfaceFluxHistory";
        surfaceFluxHistoryFile.reset(new OFstream(outputPath/fileName));
        surfaceFluxHistoryFile() << "Time(s)" << " " << "dt (s)" << " " << "uStar (m/s)" << " "
                                                                        << "TWall (K)" << " "
                                                                        << "qWall (K-m/s)" << endl;
    }
}

void Foam::surfaceShearStressFvPatchField::writeOutputFile(scalar uStar, scalar TWall, scalar qWall)
{
    if (Pstream::master)
    {
        const fvMesh& thisMesh = surfaceShearStressFvPatchField::patch().boundaryMesh().mesh();

        surfaceFluxHistoryFile() << thisMesh.time().timeName() << " " 
                                 << thisMesh.time().deltaT().value() << " " 
                                 << uStar << " " << TWall << " " << qWall << endl;
    }
}

void Foam::surfaceShearStressFvPatchField::write(Ostream& os) const
{
    fvPatchField<symmTensor>::write(os);
    os.writeKeyword("kappa") << kappa << token::END_STATEMENT << nl;
  //z0.writeEntry("z0", os);
  //os.writeKeyword("z0") << z0 << token::END_STATEMENT << nl;
    os.writeKeyword("beta_m") << beta_m << token::END_STATEMENT << nl;
    os.writeKeyword("beta_h") << beta_h << token::END_STATEMENT << nl;
    os.writeKeyword("gamma_m") << gamma_m << token::END_STATEMENT << nl;
    os.writeKeyword("gamma_h") << gamma_h << token::END_STATEMENT << nl;
    os.writeKeyword("alpha_h") << alpha_h << token::END_STATEMENT << nl;
    os.writeKeyword("B") << B << token::END_STATEMENT << nl;
    os.writeKeyword("nu") << nu << token::END_STATEMENT << nl;
    os.writeKeyword("fluxProfileRelationType") << fluxProfileRelTypeName << token::END_STATEMENT << nl;
    os.writeKeyword("averageType") << avgTypeName << token::END_STATEMENT << nl;
    os.writeKeyword("fluctuationModel") << fluctModelName << token::END_STATEMENT << nl;
    os.writeKeyword("interpolationScheme") << interpolationScheme << token::END_STATEMENT << nl;
    z0.writeEntry("z0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchSymmTensorField,
        surfaceShearStressFvPatchField
    );
}

// ************************************************************************* //
