within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function temperature_pd_der "time derivative of temperature_pd, derivative annotation of the public function temperature_pd"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Density d;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real d_der "time derivative of d";
  output Real T_der "time derivative of temperature_pd";
algorithm
  T_der := temperature_pd_state_der(p=p, d=d, state=setState_pd(p=p, d=d, phase=phase), p_der=p_der, d_der=d_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_pd is shared with the one in temperature_pd (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.temperature_pd_state_der\">temperature_pd_state_der</a>.
</p>
</html>"));
end temperature_pd_der;
