// Use trilinear interpolation using nearest cell center and cell-center
// gradient information.
DynamicList<label> bladeInfluenceCells_ = bladeInfluenceCellsGrouped[i].returnRegionCellIDs(point - vector(bladeProjectionRadius[i],bladeProjectionRadius[i],bladeProjectionRadius[i]),
                                                                                            point + vector(bladeProjectionRadius[i],bladeProjectionRadius[i],bladeProjectionRadius[i]));
//forAll(bladeInfluenceCells[i], m)
forAll(bladeInfluenceCells_, m)
{
  //vector disVector = (mesh_.C()[bladeInfluenceCells[i][m]] - point);
    vector disVector = (mesh_.C()[bladeInfluenceCells_[m]] - point);
    scalar dis = mag(disVector);
    if (dis <= bladeProjectionRadius[i])
    {
        // Compute the body force projection.
        scalar spreading = computeBladeProjectionFunction(disVector,r0,i,j,k);

        // Sum up this mesh cell's contribution to the integrated velocity, 
        // weighted by the projection function.
      //velocity += U_[bladeInfluenceCells[i][m]] * spreading * mesh_.V()[bladeInfluenceCells[i][m]];
        velocity += U_[bladeInfluenceCells_[m]] * spreading * mesh_.V()[bladeInfluenceCells_[m]];
    }
}
