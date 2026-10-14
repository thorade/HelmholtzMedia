within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function setState_dTX_der "time derivative of setState_dTX (derivative of the state record)"
  extends Modelica.Icons.Function;
  input Density d "Density";
  input Temperature T "Temperature";
  input MassFraction X[:] "Mass fractions";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  input Real X_der[:] "time derivative of X (not used, pure substance)";
  output ThermodynamicState state_der "time derivative of the state (phase is set to 0, it has no derivative)";
algorithm
  state_der := ThermodynamicState(
    phase=0,
    T=T_der,
    p=pressure_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, X=X, phase=phase), der_d=d_der, der_T=T_der),
    d=d_der,
    u=specificEnthalpy_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, X=X, phase=phase), der_d=d_der, der_T=T_der)
      - pressure_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, X=X, phase=phase), der_d=d_der, der_T=T_der)/density(state=setState_dTX(d=d, T=T, X=X, phase=phase))
      + pressure(state=setState_dTX(d=d, T=T, X=X, phase=phase))/density(state=setState_dTX(d=d, T=T, X=X, phase=phase))^2*d_der,
    h=specificEnthalpy_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, X=X, phase=phase), der_d=d_der, der_T=T_der),
    s=specificEntropy_dT_state_der(d=d, T=T, state=setState_dTX(d=d, T=T, X=X, phase=phase), d_der=d_der, T_der=T_der));
  annotation (Inline=true, Documentation(info="<html>
<p>
Derivative annotation of setState_dTX, see <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.setState_phX_der\">setState_phX_der</a>.
Written as one statement and inlined, so that the setState_dTX calls are shared with the original call (same arguments).
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_dTX_der;
