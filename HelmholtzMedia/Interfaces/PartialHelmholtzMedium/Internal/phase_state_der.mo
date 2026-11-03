within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function phase_state_der "time derivative of phase_state: zero, the phase is piecewise constant"
  extends Modelica.Icons.Function;
  input ThermodynamicState state;
  output Real phase_der "time derivative of phase_state";
algorithm
  phase_der := 0;
  annotation (Inline=true);
end phase_state_der;
