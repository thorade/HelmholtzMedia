within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function density_pT_state
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output Density d "Density";

algorithm
  d := density(state);

annotation (
  Inline=false,
  LateInline=true,
  derivative(noDerivative=state)=density_pT_der);
end density_pT_state;
