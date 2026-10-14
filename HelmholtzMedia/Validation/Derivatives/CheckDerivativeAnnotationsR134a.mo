within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsR134a "CheckDerivativeAnnotationsAllCases for R134a"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsR134a;
