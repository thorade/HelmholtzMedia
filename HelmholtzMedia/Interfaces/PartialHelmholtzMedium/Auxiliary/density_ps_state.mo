within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function density_ps_state "density_ps with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  output Density d "Density";
algorithm
  d := density(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=density_ps_der);
end density_ps_state;
