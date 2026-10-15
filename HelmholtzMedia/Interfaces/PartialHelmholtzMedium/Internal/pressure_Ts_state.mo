within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_Ts_state "pressure_Ts with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (T, s)";
  output AbsolutePressure p "Pressure";
algorithm
  p := pressure(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=pressure_Ts_state_der);
end pressure_Ts_state;
