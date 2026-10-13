within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificEntropy_pT_state "specificEntropy_pT with the thermodynamic state as additional input (carries the derivative annotation)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input ThermodynamicState state "state at (p, T)";
  output SpecificEntropy s "Specific entropy";
algorithm
  s := specificEntropy(state);
annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=specificEntropy_pT_der);
end specificEntropy_pT_state;
