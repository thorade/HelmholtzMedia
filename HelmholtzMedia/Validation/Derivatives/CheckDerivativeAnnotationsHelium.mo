within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeAnnotationsHelium "CheckDerivativeAnnotationsAllCases for Helium"
  extends CheckDerivativeAnnotationsAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Helium);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeAnnotationsHelium;
