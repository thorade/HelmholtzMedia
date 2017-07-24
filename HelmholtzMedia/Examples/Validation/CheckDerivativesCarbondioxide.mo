within HelmholtzMedia.Examples.Validation;
model CheckDerivativesCarbondioxide
  extends HelmholtzMedia.Examples.Validation.CheckDerivativesMedium(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Carbondioxide,
    p0 = 20.3e5,  p1 = 60e5,
    T0 = 300,     T1 = 500);
   annotation(experiment(StopTime= 1, Tolerance=1e-010));
end CheckDerivativesCarbondioxide;
