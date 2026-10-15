within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEntropy_pT_state_der "time derivative of specificEntropy_pT"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input ThermodynamicState state "state at (p, T)";
  input Real p_der "time derivative of p";
  input Real T_der "time derivative of T";
  output Real s_der "time derivative of s";
algorithm
  s_der := p_der*density_derT_p(state)/state.d^2 + T_der*specificHeatCapacityCp(state)/state.T;
annotation (Inline=true);
end specificEntropy_pT_state_der;
