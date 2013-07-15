within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_dT_state
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  output AbsolutePressure p "pressure";

algorithm
  p := pressure(state);

annotation (
  Inline=false,
  LateInline=true,
  inverse(d=density_pT_state(p=p, T=T, state=state)),
  derivative(noDerivative=state)=pressure_dT_der);
end pressure_dT_state;
