within HelmholtzMedia.Examples.Validation;
model Validate_HelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.Density d=838.025;
  Medium.Temperature T=500;
  HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.HelmholtzDerivs
                         f;

algorithm
    f := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.setHelmholtzDerivs(
                                   T=T, d=d, phase=1);

annotation (experiment(NumberOfIntervals=1));
end Validate_HelmholtzDerivatives;
