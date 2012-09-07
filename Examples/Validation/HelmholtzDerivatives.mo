within HelmholtzMedia.Examples.Validation;
function HelmholtzDerivatives
  // validate derivatives of Helmholtz energy (single phase state)
  // values for comparison are given in IAPWS-95 (Table 6)
  // http://iapws.org/relguide/IAPWS-95.htm

  package medium = HelmholtzMedia.HelmholtzFluids.Butane;
  constant medium.Density d=838.025;
  constant medium.Temperature T=500;
  output medium.EoS.HelmholtzDerivs f;

algorithm
    f := medium.EoS.setHelmholtzDerivsSecond( T=T, d=d, phase=1);

annotation (experiment(NumberOfIntervals=1));
end HelmholtzDerivatives;
