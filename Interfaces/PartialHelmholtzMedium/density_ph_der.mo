within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_ph_der "time derivative of density_ph"

  input AbsolutePressure p;
  input SpecificEnthalpy h;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real d_der "time derivative of density";

protected
  ThermodynamicState state=setState_phX(p=p, h=h, phase=phase);

algorithm
  d_der := p_der*density_derp_h(state=state)
         + h_der*density_derh_p(state=state);

annotation (inline=true);
end density_ph_der;
