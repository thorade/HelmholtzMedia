within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function temperature_pd_der "time derivative of temperature_pd"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input ThermodynamicState state "state at (p, d)";
  input Real p_der "time derivative of p";
  input Real d_der "time derivative of d";
  output Real T_der "time derivative of T";
algorithm
  T_der := (p_der - d_der*pressure_derd_T(state))/pressure_derT_d(state);
annotation (Inline=true);
end temperature_pd_der;
