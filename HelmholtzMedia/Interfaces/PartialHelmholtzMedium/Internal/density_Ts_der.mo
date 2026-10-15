within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_Ts_der "time derivative of density_Ts, derivative annotation of the public function density_Ts"
  extends Modelica.Icons.Function;
  input Temperature T;
  input SpecificEntropy s;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real T_der "time derivative of T";
  input Real s_der "time derivative of s";
  output Real d_der "time derivative of density_Ts";
algorithm
  d_der := density_Ts_state_der(T=T, s=s, state=setState_Ts(T=T, s=s, phase=phase), T_der=T_der, s_der=s_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_Ts is shared with the one in density_Ts (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.density_Ts_state_der\">density_Ts_state_der</a>.
</p>
</html>"));
end density_Ts_der;
