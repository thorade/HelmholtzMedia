within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_pT_state
  "returns specific enthalpy for given p and T"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output SpecificEnthalpy h "specific enthalpy";

algorithm
  h := specificEnthalpy(state);

annotation (
  Inline=false,
  LateInline=true,
  inverse(T=temperature_ph_state(p=p, h=h, state=state)),
  derivative(noDerivative=state)=specificEnthalpy_pT_der);
end specificEnthalpy_pT_state;
