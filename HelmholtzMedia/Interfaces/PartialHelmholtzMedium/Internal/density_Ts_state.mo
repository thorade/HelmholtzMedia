within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_Ts_state "density_Ts with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (T, s)";
  output Density d "Density";
algorithm
  d := density(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=density_Ts_der);
end density_Ts_state;
