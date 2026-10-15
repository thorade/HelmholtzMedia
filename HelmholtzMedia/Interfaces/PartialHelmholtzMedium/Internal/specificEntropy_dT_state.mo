within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEntropy_dT_state "specificEntropy_dT with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
  input ThermodynamicState state "state at (d, T)";
  output SpecificEntropy s "Specific entropy";
algorithm
  s := specificEntropy(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=specificEntropy_dT_state_der);
end specificEntropy_dT_state;
