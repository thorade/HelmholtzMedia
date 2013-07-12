within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_dT_der "time derivative of specificEnthalpy_dT"

  input Density d "Density";
  input Temperature T "Temperature";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real der_d "Density derivative";
  input Real der_T "Temperature derivative";
  output Real der_h "Time derivative of enthalpy";

protected
  ThermodynamicState state=setState_dT( d=d, T=T, phase=phase);

algorithm
  der_h := der_d*specificEnthalpy_derd_T(state=state)
         + der_T*specificEnthalpy_derT_d(state=state);
end specificEnthalpy_dT_der;
