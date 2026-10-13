within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function temperature_pd_state "temperature_pd with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input ThermodynamicState state "state at (p, d)";
  output Temperature T "Temperature";
algorithm
  T := temperature(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=temperature_pd_der);
end temperature_pd_state;
