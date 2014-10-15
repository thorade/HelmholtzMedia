within HelmholtzMedia.Examples.Validation;
model TestDifferentiationWater
  extends TestDifferentiationMedium(redeclare package Medium =
        Modelica.Media.Water.WaterIF97_ph);
end TestDifferentiationWater;
