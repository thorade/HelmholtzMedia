within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function phase_state "returns the phase of a state as Real (1 single phase, 2 two-phase)"
  extends Modelica.Icons.Function;
  input ThermodynamicState state;
  output Real phase "1 single phase, 2 two-phase";
algorithm
  phase := state.phase;
  annotation (Inline=false, LateInline=true, Documentation(info="<html>
<p>
Used in BaseProperties for state.phase. Like density_ph_state etc., it receives the state as a visible argument and
is not inlined, so that Dymola reuses the state solve of BaseProperties instead of solving the state again.
The Real output makes the comparison in BaseProperties a Real relation, so a change of phase is an event.
</p>
</html>"));
end phase_state;
