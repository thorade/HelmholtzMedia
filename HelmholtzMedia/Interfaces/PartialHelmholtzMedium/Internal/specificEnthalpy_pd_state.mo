within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_pd_state "specificEnthalpy_pd with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input ThermodynamicState state "state at (p, d)";
  output SpecificEnthalpy h "Specific enthalpy";
algorithm
  h := specificEnthalpy(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=specificEnthalpy_pd_der);
end specificEnthalpy_pd_state;
