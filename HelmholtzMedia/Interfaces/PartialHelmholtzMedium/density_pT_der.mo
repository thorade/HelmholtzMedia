within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_pT_der "time derivative of density_pT"

  input AbsolutePressure p;
  input Temperature T;
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real T_der "time derivative of temperature";
  output Real d_der "time derivative of density";

algorithm
  d_der := p_der*density_derp_T(state=state)
         + T_der*density_derT_p(state=state);

annotation (Inline=true);
end density_pT_der;
