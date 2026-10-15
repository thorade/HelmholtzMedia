within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEntropy_ph_der "time derivative of specificEntropy_ph, derivative annotation of the public function specificEntropy_ph"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input SpecificEnthalpy h;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real h_der "time derivative of h";
  output Real s_der "time derivative of specificEntropy_ph";
algorithm
  s_der := specificEntropy_ph_state_der(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase), p_der=p_der, h_der=h_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_ph is shared with the one in specificEntropy_ph (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEntropy_ph_state_der\">specificEntropy_ph_state_der</a>.
</p>
</html>"));
end specificEntropy_ph_der;
