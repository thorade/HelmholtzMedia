within HelmholtzMedia.Examples.Tests.Validation;
model Validate_HelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.Density d=838.025;
  Medium.Temperature T=500;
  Medium.HelmholtzDerivs f;

algorithm
    f := Medium.setHelmholtzDerivs(T=T, d=d, phase=1);

annotation (experiment(NumberOfIntervals=1));
end Validate_HelmholtzDerivatives;
