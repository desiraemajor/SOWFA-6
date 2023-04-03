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
    Foam::cartesianMeshSearch

Description
    Divides a region of the mesh into a regular Cartesian array of subvolumes
    and figures out which mesh cell centers lie within each subvolume.  This
    can be used to speed up the search for a grid cell.  For example, this is
    used to find sampling and body force projection grid cells in the actuator
    line model.

SourceFiles
    cartesianMeshSearch.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshSearch_H
#define cartesianMeshSearch_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshSearch
{
    // Private data

        //- Reference to mesh
      //const fvMesh& mesh_;


        //- i,j,k list of cells grouped into Cartesian subvolumes.
        List<List<List<DynamicList<label> > > > cellIDListGrouped_;

        //- i,j,k list of bounding boxes of the Cartesian subvolumes.
        List<List<List<boundBox> >  > boundBoxListGrouped_;

        //- Number of subvolumes in each direction.
        int ni_;
        int nj_;
        int nk_;

        //- Dimensions of subvolumes
        vector subVolumeSize_;

        //- Minimum and maximum point on the Cartesian block.
        vector blockMin_, blockMax_;



    // Private Member Functions

        //- Find cells per plane 
        void findPlanes_();


public:

    //- Declare name of the class and its debug switch
    ClassName("cartesianMeshSearch");



    // Constructors

        //- Default constructor.
        cartesianMeshSearch();

        //- Fully initializing constructor.
        /*
        cartesianMeshSearch
        (
            const fvMesh& mesh,
            const vector blockMin,
            const vector blockMax,
            const vector subBlockSize,
            const DynamicList<label>& cellIDList
        );
        */


    // Destructor
    virtual ~cartesianMeshSearch();


    // Member functions

        //- Access

            //- Do the search and organization of grid cells.
            void update
            (
                const fvMesh& mesh,
                const DynamicList<label>& cellIDList,
                const vector blockMin,
                const vector blockMax,
                vector subVolumeSize
            );

            
            //- Return the cell IDs for the mesh in a sub region.
            DynamicList<label> returnRegionCellIDs
            (
                const vector regionMin, 
                const vector regionMax
            );
            

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
