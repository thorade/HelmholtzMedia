within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetStateHelium "CheckDerivativeSetStateAllCases for Helium"
  extends CheckDerivativeSetStateAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Helium);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeSetStateHelium;
