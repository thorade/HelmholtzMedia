within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsButane "CheckDerivativeAnnotationsAllCases for Butane"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsButane;
