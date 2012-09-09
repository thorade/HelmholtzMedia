within HelmholtzMedia.Examples.Validation;
model HelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package medium = HelmholtzMedia.HelmholtzFluids.Butane;
  parameter medium.Density d=838.025;
  parameter medium.Temperature T=500;
  medium.EoS.HelmholtzDerivs f=medium.EoS.setHelmholtzDerivsSecond(T=T, d=d, phase=1);

annotation (experiment(NumberOfIntervals=1));
end HelmholtzDerivatives;
