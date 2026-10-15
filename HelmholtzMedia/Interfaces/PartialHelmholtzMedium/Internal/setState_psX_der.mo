within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function setState_psX_der "time derivative of setState_psX (derivative of the state record)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input MassFraction X[:] "Mass fractions";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real s_der "time derivative of s";
  input Real X_der[:] "time derivative of X (not used, pure substance)";
  output ThermodynamicState state_der "time derivative of the state (phase is set to 0, it has no derivative)";
protected
  ThermodynamicState state "the state, evaluated once";
  HelmholtzCache cache "evaluated once";
  Real d_der;
  Real h_der;
  Real T_der;
algorithm
  if (p_der == 0 and s_der == 0) then
    // zero seed (e.g. in the analytic Jacobian of a nonlinear system, for a state that does not depend on
    // the iteration variables): skip the evaluation of the state
    state_der := ThermodynamicState(phase=0, d=0, T=0, p=0, h=0, u=0, s=0);
    return;
  end if;
  state := setState_psX(p=p, s=s, X=X, phase=phase);
  cache := setHelmholtzCache(d=state.d, T=state.T, phase=state.phase);
  T_der := if state.phase == 2 then p_der*jouleThomsonCoefficient_cache(d=state.d, T=state.T, phase=state.phase, cache=cache)
    else (-p_der*state.T*density_derT_p_cache(d=state.d, T=state.T, phase=state.phase, cache=cache)/state.d^2 + s_der*state.T)/specificHeatCapacityCp_cache(d=state.d, T=state.T, phase=state.phase, cache=cache);
  d_der := p_der/velocityOfSound_cache(d=state.d, T=state.T, phase=state.phase, cache=cache)^2 + s_der*state.T*density_derh_p_cache(d=state.d, T=state.T, phase=state.phase, cache=cache);
  h_der := p_der/state.d + s_der*state.T;
  state_der := ThermodynamicState(
    phase=0,
    T=T_der,
    p=p_der,
    d=d_der,
    u=h_der - p_der/state.d + p/state.d^2*d_der,
    h=h_der,
    s=s_der);
  annotation (Documentation(info="<html>
<p>
Derivative annotation of setState_psX, see <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.setState_phX_der\">setState_phX_der</a>.
Not inlined: the state and the HelmholtzCache are evaluated once (not at all if both derivatives of the inputs are zero), and the partial derivatives are computed by the *_cache workers.
Dymola does not inline derivative functions inside the analytic Jacobians of nonlinear systems, where an inlined version evaluated the state several times per call.
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_psX_der;
