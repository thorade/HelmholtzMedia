within HelmholtzFluids.Examples;
model State_dT_parameter "calculate state record from dT input"

  package medium = HelmholtzFluids.R134a;

  parameter medium.Density d=1000;
  parameter medium.Temperature T=400;

  medium.ThermodynamicState inletState;

  // medium.MassFraction x;
  // medium.SurfaceTension sigma;
  medium.DynamicViscosity eta;

equation
  inletState=medium.setState_dTX(d=d, T=T, phase=0, X={1});

  // x=medium.vapourQuality(inletState);
  // sigma=medium.surfaceTension(medium.setSat_T(T=T));
  eta=medium.dynamicViscosity(inletState);

end State_dT_parameter;
