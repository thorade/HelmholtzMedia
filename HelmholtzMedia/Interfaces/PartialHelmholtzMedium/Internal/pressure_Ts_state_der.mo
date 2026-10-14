within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_Ts_state_der "time derivative of pressure_Ts"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (T, s)";
  input Real T_der "time derivative of T";
  input Real s_der "time derivative of s";
  output Real p_der "time derivative of p";
algorithm
  p_der := if state.phase == 2 then T_der*pressure_derT_d(state)
    else (-T_der*specificHeatCapacityCp(state)/state.T + s_der)*state.d^2/density_derT_p(state);
annotation (Inline=true);
end pressure_Ts_state_der;
