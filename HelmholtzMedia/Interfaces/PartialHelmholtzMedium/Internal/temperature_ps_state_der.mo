within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function temperature_ps_state_der "time derivative of temperature_ps"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  input Real p_der "time derivative of p";
  input Real s_der "time derivative of s";
  output Real T_der "time derivative of T";
algorithm
  T_der := if state.phase == 2 then p_der*jouleThomsonCoefficient(state)
    else (-p_der*state.T*density_derT_p(state)/state.d^2 + s_der*state.T)/specificHeatCapacityCp(state);
annotation (Inline=true);
end temperature_ps_state_der;
