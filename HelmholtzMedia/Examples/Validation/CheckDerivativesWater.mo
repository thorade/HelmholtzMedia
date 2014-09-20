within HelmholtzMedia.Examples.Validation;
model CheckDerivativesWater
  extends CheckDerivativesMedium(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p0 = 20.3e5,  p1 = 60e5,
    T0 = 800,     T1 = 1000);
  annotation (experiment(StopTime=1, Tolerance=1e-010));
end CheckDerivativesWater;
