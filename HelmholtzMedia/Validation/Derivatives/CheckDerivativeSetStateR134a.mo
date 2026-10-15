within HelmholtzMedia.Validation.Derivatives;
model CheckDerivativeSetStateR134a "CheckDerivativeSetStateAllCases for R134a"
  extends CheckDerivativeSetStateAllCases(redeclare package Medium = HelmholtzMedia.HelmholtzFluids.R134a);
  annotation (experiment(StopTime=1, Tolerance=1e-10));
end CheckDerivativeSetStateR134a;
