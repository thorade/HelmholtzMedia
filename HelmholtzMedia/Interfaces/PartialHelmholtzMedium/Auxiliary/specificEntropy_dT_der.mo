within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificEntropy_dT_der "time derivative of specificEntropy_dT"
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
  input ThermodynamicState state "state at (d, T)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real s_der "time derivative of s";
algorithm
  s_der := -d_der*pressure_derT_d(state)/state.d^2 + T_der*specificHeatCapacityCv(state)/state.T;
annotation (Inline=true);
end specificEntropy_dT_der;
