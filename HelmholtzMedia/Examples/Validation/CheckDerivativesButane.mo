within HelmholtzMedia.Examples.Validation;
model CheckDerivativesButane
  extends HelmholtzMedia.Examples.Validation.CheckDerivativesMedium(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Butane,
    p0 = 20.3e5,  p1 = 60e5,
    T0 = 300,     T1 = 500);
   annotation(experiment(StopTime= 1, Tolerance=1e-010));
end CheckDerivativesButane;
