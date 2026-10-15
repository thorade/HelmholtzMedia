within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetStateAmmonia "CheckDerivativeSetStateAllCases for Ammonia"
  extends CheckDerivativeSetStateAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeSetStateAmmonia;
