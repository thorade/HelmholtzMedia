within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_ps_state "specificEnthalpy_ps with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  output SpecificEnthalpy h "Specific enthalpy";
algorithm
  h := specificEnthalpy(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=specificEnthalpy_ps_state_der);
end specificEnthalpy_ps_state;
