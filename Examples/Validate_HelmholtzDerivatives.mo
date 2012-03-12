within HelmholtzMedia.Examples;
model Validate_HelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package medium = HelmholtzMedia.HelmholtzFluids.Butane;
  medium.Density d=838.025;
  medium.Temperature T=500;
  medium.HelmholtzDerivs f;

algorithm
    f := medium.setHelmholtzDerivs(T=T, d=d, phase=1);
  annotation (experiment(StopTime=10, NumberOfIntervals=1),
                                       __Dymola_experimentSetupOutput);
end Validate_HelmholtzDerivatives;
