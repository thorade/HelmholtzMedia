within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_dT_der "time derivative of specificEnthalpy_dT, derivative annotation of the public function specificEnthalpy_dT"
  extends Modelica.Icons.Function;
  input Density d;
  input Temperature T;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real der_h "time derivative of specificEnthalpy_dT";
algorithm
  der_h := specificEnthalpy_dT_state_der(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase), der_d=d_der, der_T=T_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_dT is shared with the one in specificEnthalpy_dT (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEnthalpy_dT_state_der\">specificEnthalpy_dT_state_der</a>.
</p>
</html>"));
end specificEnthalpy_dT_der;
