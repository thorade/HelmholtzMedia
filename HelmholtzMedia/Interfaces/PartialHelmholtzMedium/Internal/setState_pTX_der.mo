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
protected
  ThermodynamicState state=setState_pTX(p=p, T=T, X=X, phase=phase) "the state, evaluated once";
  HelmholtzCache cache=setHelmholtzCache(d=state.d, T=state.T, phase=state.phase) "evaluated once";
  Real d_der;
  Real h_der;
  Real ddTp;
  Real cp;
algorithm
  ddTp := density_derT_p_cache(d=state.d, T=state.T, phase=state.phase, cache=cache);
  cp := specificHeatCapacityCp_cache(d=state.d, T=state.T, phase=state.phase, cache=cache);
  d_der := p_der*density_derp_T_cache(d=state.d, T=state.T, phase=state.phase, cache=cache) + T_der*ddTp;
  h_der := p_der*isothermalThrottlingCoefficient_cache(d=state.d, T=state.T, phase=state.phase, cache=cache) + T_der*cp;
  state_der := ThermodynamicState(
    phase=0,
    T=T_der,
    p=p_der,
    d=d_der,
    u=h_der - p_der/state.d + state.p/state.d^2*d_der,
    h=h_der,
    s=p_der*ddTp/state.d^2 + T_der*cp/state.T);
  annotation (Documentation(info="<html>
<p>
Derivative annotation of setState_pTX, see <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.setState_phX_der\">setState_phX_der</a>.
Not inlined: the state and the HelmholtzCache are evaluated once, and the partial derivatives are computed by the *_cache workers.
Dymola does not inline derivative functions inside the analytic Jacobians of nonlinear systems, where an inlined version evaluated the state several times per call.
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_pTX_der;
