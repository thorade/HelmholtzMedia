within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_ps_der "time derivative of density_ps"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  input Real p_der "time derivative of p";
  input Real s_der "time derivative of s";
  output Real d_der "time derivative of d";
algorithm
  d_der := p_der/velocityOfSound(state)^2 + s_der*state.T*density_derh_p(state);
annotation (Inline=true);
end density_ps_der;
