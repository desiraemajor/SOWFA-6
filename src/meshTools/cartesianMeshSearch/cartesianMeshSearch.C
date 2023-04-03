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

#include "cartesianMeshSearch.H"



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Default constructor.
Foam::cartesianMeshSearch::cartesianMeshSearch()
:
ni_(1),
nj_(1),
nk_(1),
subVolumeSize_(vector(1.0, 1.0, 1.0)),
blockMin_(vector(0.0, 0.0, 0.0)),
blockMax_(vector(1.0, 1.0, 1.0))
{
    // Do nothing.
}


// Constructor that fully initializes the Cartesian mesh search.
/*
Foam::cartesianMeshSearch::cartesianMeshSearch
(
    const fvMesh& mesh,
    const vector blockMin,
    const vector blockMax,
    const vector subBlockSize,
    const DynamicList<label>& cellIDList
)
:
    mesh_(mesh),
    blockMin_(blockMin),
    blockMax_(blockMax),
    subBlockSize_(subBlockSize),
    initialized_(false)
{
    // Create the Cartesian cell grouping.
    update
    (
        mesh_,
        cellIDList,
        blockMin_,
        blockMax_,
        subBlockSize_
    );

    initialized_ = true;
}
*/


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cartesianMeshSearch::~cartesianMeshSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cartesianMeshSearch::update
(
    const fvMesh& mesh,
    const DynamicList<label>& cellIDList,
    const vector blockMin,
    const vector blockMax,
    vector subVolumeSize
)
{
    // Set the internal private variables.
    subVolumeSize_ = subVolumeSize;
    blockMin_ = blockMin;
    blockMax_ = blockMax;

    // Figure out the subblocks.
    ni_ = round((blockMax_.x() - blockMin_.x()) / subVolumeSize_.x());
    nj_ = round((blockMax_.y() - blockMin_.y()) / subVolumeSize_.y());
    nk_ = round((blockMax_.z() - blockMin_.z()) / subVolumeSize_.z());

    subVolumeSize_.x() = (blockMax_.x() - blockMin_.x()) / ni_;
    subVolumeSize_.y() = (blockMax_.y() - blockMin_.y()) / nj_;
    subVolumeSize_.z() = (blockMax_.z() - blockMin_.z()) / nk_;


    // Layout the grouped boundBox list.
    {
        List<boundBox> kList_;
        kList_.setSize(nk_);
        List<List<boundBox> > jkList_;
        jkList_.setSize(nj_);
        List<List<List<boundBox> > > ijkList_;
        ijkList_.setSize(ni_);
        for (int j = 0; j < nj_; j++)
        {
            jkList_[j] = kList_;
        }
        for (int i = 0; i < ni_; i++)
        {
            ijkList_[i] = jkList_;
        }
        boundBoxListGrouped_ = ijkList_;
    }

    // Layout the grouped cellID list.
    {
        List<DynamicList<label> > kList_;
        kList_.setSize(nk_);
        List<List<DynamicList<label> > > jkList_;
        jkList_.setSize(nj_);
        List<List<List<DynamicList<label> > > > ijkList_;
        ijkList_.setSize(ni_);
        for (int j = 0; j < nj_; j++)
        {
            jkList_[j] = kList_;
        }
        for (int i = 0; i < ni_; i++)
        {
            ijkList_[i] = jkList_;
        }
        cellIDListGrouped_ = ijkList_;
    }

    // Load the bounding box list.
    for (int i = 0; i < ni_; i++)
    {
        for (int j = 0; j < nj_; j++)
        {
            for (int k = 0; k < nk_; k++)
            {
                vector bbMin = vector
                               (
                                   blockMin_.x() + i*subVolumeSize_.x(), 
                                   blockMin_.y() + j*subVolumeSize_.y(), 
                                   blockMin_.z() + k*subVolumeSize_.z()
                               );
                vector bbMax = vector
                               (
                                   blockMin_.x() + (i+1)*subVolumeSize_.x(), 
                                   blockMin_.y() + (j+1)*subVolumeSize_.y(), 
                                   blockMin_.z() + (k+1)*subVolumeSize_.z()
                               );
                boundBoxListGrouped_[i][j][k] = boundBox(bbMin,bbMax);

                forAll(mesh.C(),m)
                {
                    vector meshPoint = mesh.C()[m];
                    if (boundBoxListGrouped_[i][j][k].contains(meshPoint))
                    {
                        cellIDListGrouped_[i][j][k].append(m);
                    }
                }
            }
        }
    }

/*
    // Populate the grouped cellID list.
    for (int i = 0; i < ni_; i++)
    {
        for (int j = 0; j < nj_; j++)
        {
            for (int k = 0; k < nk_; k++)
            {
              //forAll(cellIDList,m)
                forAll(mesh.C(),m)
                {
                  //vector meshPoint = mesh.C()[cellIDList[m]];
                    vector meshPoint = mesh.C()[m];
                  //Info << meshPoint << endl;
                  //Info << boundBoxListGrouped_[i][j][k] << tab << meshPoint << endl;
                    if (boundBoxListGrouped_[i][j][k].contains(meshPoint))
                    {
                      //cellIDListGrouped_[i][j][k].append(cellIDList[m]);
                        cellIDListGrouped_[i][j][k].append(m);
                    }
                }
            }
        }
    }
    */

    /*
    Info << "In update..." << endl;
    Info << "blockMin = " << blockMin_ << endl;
    Info << "blockMax = " << blockMax_ << endl;
    Info << "subVolumeSize = " << subVolumeSize_ << endl;
    Info << "ni, nj, nk = " << ni_ << "," << nj_ << "," << nk_ << endl;
    Info << "boundBoxListGrouped = " << boundBoxListGrouped_ << endl;
    Info << "cellIDList.size() = " << cellIDList.size() << endl;
    Info << "cellIDListGrouped.size() = " << cellIDListGrouped_.size() << endl;
    */
    
}


DynamicList<label> cartesianMeshSearch::returnRegionCellIDs
(
    const vector regionMin, 
    const vector regionMax
)
{
  //Info << "In returnRegionCellIDs..." << endl;
   
    label iMin = max(floor((regionMin.x() - blockMin_.x()) / subVolumeSize_.x()), 0);
    label iMax = min(floor((regionMax.x() - blockMin_.x()) / subVolumeSize_.x()), ni_-1);

    label jMin = max(floor((regionMin.y() - blockMin_.y()) / subVolumeSize_.y()), 0);
    label jMax = min(floor((regionMax.y() - blockMin_.y()) / subVolumeSize_.y()), nj_-1);

    label kMin = max(floor((regionMin.z() - blockMin_.z()) / subVolumeSize_.z()), 0);
    label kMax = min(floor((regionMax.z() - blockMin_.z()) / subVolumeSize_.z()), nk_-1);

    /*
    Info << "iMin:iMax = " << iMin << ":" << iMax << endl;
    Info << "jMin:jMax = " << jMin << ":" << jMax << endl;
    Info << "kMin:kMax = " << kMin << ":" << kMax << endl;

    Info << "xMin:xMax = " << blockMin_.x() + iMin*subVolumeSize_.x() << ":" << blockMin_.x() + (iMax+1)*subVolumeSize_.x() << endl;
    Info << "yMin:yMax = " << blockMin_.y() + jMin*subVolumeSize_.y() << ":" << blockMin_.y() + (jMax+1)*subVolumeSize_.y() << endl;
    Info << "zMin:zMax = " << blockMin_.z() + kMin*subVolumeSize_.z() << ":" << blockMin_.z() + (kMax+1)*subVolumeSize_.z() << endl;

    Info << "regionMin = " << regionMin << endl;
    Info << "regionMax = " << regionMax << endl;
    */

    DynamicList<label> regionCellIDList;
    for (int i = iMin; i <= iMax; i++)
    {
        for (int j = jMin; j <= jMax; j++)
        {
            for (int k = kMin; k <= kMax; k++)
            {
              //regionCellIDList.append(cellIDListGrouped_[i][j][k]);
                forAll(cellIDListGrouped_[i][j][k],m)
                {
                    regionCellIDList.append(cellIDListGrouped_[i][j][k][m]);
                }
            }
        }
    }

  //Info << "regionCellIDList.size() = " << regionCellIDList.size() << endl;
  //Info << "regionCellIDList = " << regionCellIDList << endl;
    
    return regionCellIDList;
}


// ************************************************************************* //
