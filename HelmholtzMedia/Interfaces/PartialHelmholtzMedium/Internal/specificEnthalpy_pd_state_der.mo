within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_pd_state_der "time derivative of specificEnthalpy_pd"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input ThermodynamicState state "state at (p, d)";
  input Real p_der "time derivative of p";
  input Real d_der "time derivative of d";
  output Real h_der "time derivative of h";
algorithm
  h_der := d_der*specificEnthalpy_derd_T(state)
    + (p_der - d_der*pressure_derd_T(state))/pressure_derT_d(state)*specificEnthalpy_derT_d(state);
annotation (Inline=true);
end specificEnthalpy_pd_state_der;
