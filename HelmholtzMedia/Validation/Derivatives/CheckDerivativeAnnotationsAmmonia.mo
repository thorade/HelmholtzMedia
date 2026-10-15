within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsAmmonia "CheckDerivativeAnnotationsAllCases for Ammonia"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsAmmonia;
