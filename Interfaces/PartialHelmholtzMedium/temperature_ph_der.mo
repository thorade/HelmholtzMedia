within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_ph_der "time derivative of temperature_ph"

  input AbsolutePressure p;
  input SpecificEnthalpy h;
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real T_der "time derivative of temperature";

algorithm
  T_der := p_der*temperature_derp_h(state=state)
         + h_der*temperature_derh_p(state=state);

annotation (Inline=true);
end temperature_ph_der;
