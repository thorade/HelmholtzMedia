within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_idd "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_ideal_delta_delta
    "ideal part of dimensionless Helmholtz energy";

algorithm
    f_ideal_delta_delta := -1.0/(delta*delta);
end f_idd;
