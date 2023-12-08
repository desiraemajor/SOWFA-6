    // Solve the momentum equation
    
    // Set up the solution matrix
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      + fvc::div(Rwall)
      - Coriolis.force()
      - momentumGeoMesoTerm.force()
      - momentumSpongeLayers.force()
      - momentumPerturbationZones.source()
     ==
        (1.0 - scalar(fvOptionSourceStabilization)) * fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p_rgh)
                  + Boussinesq.buoyancyTerm()
                  + scalar(fvOptionSourceStabilization) * ((fvc::flux(fvOption(U))/mesh.magSf()) & mesh.magSf())
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

    if (pimple.finalPimpleIter() && pimple.finalPISOIter())
    {
        momentumPerturbationZones.perturb(U);
    }
