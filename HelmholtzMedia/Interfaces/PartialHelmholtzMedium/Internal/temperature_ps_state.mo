within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function temperature_ps_state "temperature_ps with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  output Temperature T "Temperature";
algorithm
  T := temperature(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=temperature_ps_state_der);
end temperature_ps_state;
