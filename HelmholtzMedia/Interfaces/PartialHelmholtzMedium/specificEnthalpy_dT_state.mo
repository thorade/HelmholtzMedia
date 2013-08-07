within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_dT_state
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output SpecificEnthalpy h "SpecificEnthalpy";

algorithm
  h := specificEnthalpy(state);

annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=specificEnthalpy_dT_der);
end specificEnthalpy_dT_state;
