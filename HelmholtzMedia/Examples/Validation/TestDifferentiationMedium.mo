within HelmholtzMedia.Examples.Validation;
model TestDifferentiationMedium
  // This example was taken from https://trac.modelica.org/Modelica/ticket/1575
  replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium;
  parameter Modelica.SIunits.Volume V=1;

  Medium.ThermodynamicState state;
  Modelica.SIunits.Density d;
  Modelica.SIunits.Mass m;
  Modelica.SIunits.Pressure p(stateSelect=StateSelect.always);
  Modelica.SIunits.SpecificInternalEnergy u;
  Modelica.SIunits.SpecificEnthalpy h(stateSelect=StateSelect.always);

equation
  m = V*d;
  state = Medium.setState_phX(p,h,Medium.reference_X);
  u = Medium.specificInternalEnergy(state);
  d = Medium.density(state);
  der(m) = time;
  der(m*u) = time;
end TestDifferentiationMedium;
