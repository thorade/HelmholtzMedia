within HelmholtzMedia.Examples.Validation;
model TestDifferentiationButane
  extends TestDifferentiationMedium(redeclare package Medium =
        HelmholtzMedia.HelmholtzFluids.Butane);
end TestDifferentiationButane;
