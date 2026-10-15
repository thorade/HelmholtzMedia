within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEntropy_pT_der "time derivative of specificEntropy_pT, derivative annotation of the public function specificEntropy_pT"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Temperature T;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real T_der "time derivative of T";
  output Real s_der "time derivative of specificEntropy_pT";
algorithm
  s_der := specificEntropy_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, phase=phase), p_der=p_der, T_der=T_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_pTX is shared with the one in specificEntropy_pT (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEntropy_pT_state_der\">specificEntropy_pT_state_der</a>.
</p>
</html>"));
end specificEntropy_pT_der;
