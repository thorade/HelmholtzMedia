within HelmholtzFluids.Examples;
model State_dT_parameter "calculate state record from dT input"

  package medium = HelmholtzFluids.R134a;

  parameter medium.Density d=228;
  parameter medium.Temperature T=273.15;

  medium.ThermodynamicState inlet;

  // medium.MassFraction x;
  // medium.SurfaceTension sigma;

equation
  inlet=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // x=medium.vapourQuality(inlet);
  // sigma=medium.surfaceTension(medium.setSat_T(T=T));

end State_dT_parameter;
