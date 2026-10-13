within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_Ts_der "time derivative of density_Ts"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (T, s)";
  input Real T_der "time derivative of T";
  input Real s_der "time derivative of s";
  output Real d_der "time derivative of d";
algorithm
  d_der := (T_der*specificHeatCapacityCv(state)/state.T - s_der)*state.d^2/pressure_derT_d(state);
annotation (Inline=true);
end density_Ts_der;
