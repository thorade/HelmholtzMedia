within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_pT_der "time derivative of specificEnthalpy_pT, derivative annotation of the public function specificEnthalpy_pT"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Temperature T;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real T_der "time derivative of T";
  output Real h_der "time derivative of specificEnthalpy_pT";
algorithm
  h_der := specificEnthalpy_pT_state_der(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase), p_der=p_der, T_der=T_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_pT is shared with the one in specificEnthalpy_pT (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEnthalpy_pT_state_der\">specificEnthalpy_pT_state_der</a>.
</p>
</html>"));
end specificEnthalpy_pT_der;
