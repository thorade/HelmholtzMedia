within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function setState_pTX_der "time derivative of setState_pTX (derivative of the state record)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Temperature T "Temperature";
  input MassFraction X[:] "Mass fractions";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real T_der "time derivative of T";
  input Real X_der[:] "time derivative of X (not used, pure substance)";
  output ThermodynamicState state_der "time derivative of the state (phase is set to 0, it has no derivative)";
algorithm
  state_der := ThermodynamicState(
    phase=0,
    T=T_der,
    p=p_der,
    d=density_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, X=X, phase=phase), p_der=p_der, T_der=T_der),
    u=specificEnthalpy_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, X=X, phase=phase), p_der=p_der, T_der=T_der)
      - p_der/density(state=setState_pTX(p=p, T=T, X=X, phase=phase))
      + pressure(state=setState_pTX(p=p, T=T, X=X, phase=phase))/density(state=setState_pTX(p=p, T=T, X=X, phase=phase))^2*density_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, X=X, phase=phase), p_der=p_der, T_der=T_der),
    h=specificEnthalpy_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, X=X, phase=phase), p_der=p_der, T_der=T_der),
    s=specificEntropy_pT_state_der(p=p, T=T, state=setState_pTX(p=p, T=T, X=X, phase=phase), p_der=p_der, T_der=T_der));
  annotation (Inline=true, Documentation(info="<html>
<p>
Derivative annotation of setState_pTX, see <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.setState_phX_der\">setState_phX_der</a>.
Written as one statement and inlined, so that the setState_pTX calls are shared with the original call (same arguments).
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_pTX_der;
