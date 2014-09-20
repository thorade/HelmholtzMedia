within HelmholtzMedia.Examples.Validation;
model CheckDerivativesHelium
  extends HelmholtzMedia.Examples.Validation.CheckDerivativesMedium(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Helium,
    p0 = 2.3e5,  p1 = 6e5,
    T0 = 3.5,    T1 = 6);
   annotation(experiment(StopTime= 1, Tolerance=1e-010));
end CheckDerivativesHelium;
