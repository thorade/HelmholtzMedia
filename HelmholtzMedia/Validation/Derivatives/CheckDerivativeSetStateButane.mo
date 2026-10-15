within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetStateButane "CheckDerivativeSetStateAllCases for Butane"
  extends CheckDerivativeSetStateAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeSetStateButane;
