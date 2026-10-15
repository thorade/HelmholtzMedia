within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsR32 "CheckDerivativeAnnotationsAllCases for R32"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R32);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsR32;
