within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEntropy_dT_der "time derivative of specificEntropy_dT, derivative annotation of the public function specificEntropy_dT"
  extends Modelica.Icons.Function;
  input Density d;
  input Temperature T;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real s_der "time derivative of specificEntropy_dT";
algorithm
  s_der := specificEntropy_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, phase=phase), d_der=d_der, T_der=T_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_dTX is shared with the one in specificEntropy_dT (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEntropy_dT_state_der\">specificEntropy_dT_state_der</a>.
</p>
</html>"));
end specificEntropy_dT_der;
