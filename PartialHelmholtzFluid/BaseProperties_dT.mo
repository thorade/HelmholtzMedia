within HelmholtzFluids.PartialHelmholtzFluid;
model BaseProperties_dT
  "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium"
 SpecificEntropy s;

equation
  R = Modelica.Constants.R/fluidConstants[1].molarMass;
  MM = fluidConstants[1].molarMass;
algorithm
  state :=setState_dTX(d=d, T=T);
  sat:=setSat_T(T=T);
  p :=state.p;
  h :=state.h;
  s :=state.s;
equation
  u = h - p/d;
end BaseProperties_dT;
