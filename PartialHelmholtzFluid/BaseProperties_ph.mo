within HelmholtzFluids.PartialHelmholtzFluid;
model BaseProperties_ph
  "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium"
 SpecificEntropy s;

equation
  R = Modelica.Constants.R/fluidConstants[1].molarMass;
  MM = fluidConstants[1].molarMass;
algorithm
  state :=setState_phX(p=p, h=h);
  T :=state.T;
  sat:=setSat_T(T=T);
  d :=state.d;
  s :=state.s;
equation
  u = h - p/d;
end BaseProperties_ph;
