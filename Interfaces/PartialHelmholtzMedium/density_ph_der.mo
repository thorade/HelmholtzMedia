within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_ph_der "time derivative of density_ph"

  input AbsolutePressure p;
  input SpecificEnthalpy h;
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real d_der "time derivative of density";

algorithm
  d_der := p_der*density_derp_h(state=state)
         + h_der*density_derh_p(state=state);

annotation (Inline=true);
end density_ph_der;
