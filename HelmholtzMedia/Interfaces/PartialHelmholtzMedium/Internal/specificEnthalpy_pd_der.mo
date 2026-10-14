within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_pd_der "time derivative of specificEnthalpy_pd, derivative annotation of the public function specificEnthalpy_pd"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Density d;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real d_der "time derivative of d";
  output Real h_der "time derivative of specificEnthalpy_pd";
algorithm
  h_der := specificEnthalpy_pd_state_der(p=p, d=d, state=setState_pd(p=p, d=d, phase=phase), p_der=p_der, d_der=d_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_pd is shared with the one in specificEnthalpy_pd (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEnthalpy_pd_state_der\">specificEnthalpy_pd_state_der</a>.
</p>
</html>"));
end specificEnthalpy_pd_der;
