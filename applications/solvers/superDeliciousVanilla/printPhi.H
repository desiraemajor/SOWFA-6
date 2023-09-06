Info << "phi info: " << endl;
forAll(phiPatches,j)
{
    scalar phiPatchGlobal = 0.0;
    if (Pstream::myProcNo() == phiFaceProc[j])
    {
        phiPatchGlobal = phi.boundaryField()[phiPatchID[j]][phiFaceID[j]];
    }
    reduce(phiPatchGlobal,sumOp<scalar>());
    
    Info << phiPatches[j] << ":" << endl;
    Info << "   " << phiPatchGlobal << endl;
}
