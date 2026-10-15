within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_pT_der "time derivative of density_pT, derivative annotation of the public function density_pT"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Temperature T;
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real T_der "time derivative of T";
  output Real d_der "time derivative of density_pT";
algorithm
  d_der := density_pT_state_der(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase), p_der=p_der, T_der=T_der);
  annotation (Inline=true, Documentation(info="<html>
<p>
Inlined, so that the call of setState_pT is shared with the one in density_pT (same scalar arguments).
The derivative itself is computed in <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.density_pT_state_der\">density_pT_state_der</a>.
</p>
</html>"));
end density_pT_der;
