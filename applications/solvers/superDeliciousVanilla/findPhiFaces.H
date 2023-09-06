List<point> phiFace(2,point(Zero));

phiFace[0].x() = 5040.0;
phiFace[1].x() = 0.0;

phiFace[0].y() = 25.0;
phiFace[1].y() = 25.0;

phiFace[0].z() = 15.0;
phiFace[1].z() = 15.0;

List<point> phiFaceLoc(2,point(Zero));


List<label> phiFaceID(2,label(-1));
List<label> phiFaceProc(2,label(-1));


List<label> phiPatchID(2,label(0));
List<word> phiPatches(2);
phiPatches[0] = "west.1";
phiPatches[1] = "west.2";

List<scalar> minDist(2,1.0E9);


// Loop over all boundary patches
forAll(phi.boundaryField(), patchi)
{
    word patchName = phi.mesh().boundary()[patchi].name();
    forAll(phiPatches,j)
    {
        // if this patch has one of the specified patch names, then find 
        // the patch face at the specified location.
        if (patchName == phiPatches[j])
        {
            phiPatchID[j] = patchi;
            scalar minDistLocal = 1.0E9;
            forAll(phi.boundaryField()[patchi],i)
            {
                vector faceLoc = mesh.boundary()[patchi].Cf()[i];
                scalar dist = Foam::sqrt(Foam::sqr(phiFace[j].y() - faceLoc.y()) + Foam::sqr(phiFace[j].z() - faceLoc.z()));

                if (dist < minDistLocal)
                {
                    minDistLocal = dist;
                    phiFaceID[j] = i;
                }
            }
            scalar minDistGlobal = minDistLocal;
            reduce(minDistGlobal, minOp<scalar>());
            minDist[j] = minDistGlobal;
            
            label phiFaceProcGlobal = 0;
            label phiFaceIDGlobal = 0;
            point phiFaceLocGlobal(Zero);
            if (minDistGlobal == minDistLocal)
            {
                phiFaceProcGlobal = Pstream::myProcNo();
                phiFaceIDGlobal = phiFaceID[j];
                phiFaceLocGlobal = mesh.boundary()[patchi].Cf()[phiFaceIDGlobal];
            }
            reduce(phiFaceLocGlobal, sumOp<point>());
            reduce(phiFaceProcGlobal, sumOp<label>());
            reduce(phiFaceIDGlobal, sumOp<label>());
            phiFaceLoc[j] = phiFaceLocGlobal;
            phiFaceProc[j] = phiFaceProcGlobal;
            phiFaceID[j] = phiFaceIDGlobal;
        }
    }
}


forAll(phiPatches,j)
{
    Info << "phi tracking on patch: " << phiPatchID[j] << ", " << phiPatches[j] << endl;
    Info << "face: " << phiFaceID[j] << " at location " << phiFaceLoc[j] << endl;
    Info << "on processor " << phiFaceProc[j] << endl;
}
