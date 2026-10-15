within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_ps_der "time derivative of density_ps, derivative annotation of the public function density_ps"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input SpecificEntropy s;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real s_der "time derivative of s";
  output Real d_der "time derivative of density_ps";
algorithm
  d_der := density_ps_state_der(p=p, s=s, state=setState_ps(p=p, s=s, phase=phase), p_der=p_der, s_der=s_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_ps is shared with the one in density_ps (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.density_ps_state_der\">density_ps_state_der</a>.
</p>
</html>"));
end density_ps_der;
