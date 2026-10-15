within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function setState_phX_der "time derivative of setState_phX (derivative record)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific enthalpy";
  input MassFraction X[:] "Mass fractions";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input Real p_der "time derivative of p";
  input Real h_der "time derivative of h";
  input Real X_der[:] "time derivative of X (not used, pure substance)";
  output ThermodynamicState state_der "time derivative of the state (phase is set to 0, it has no derivative)";
protected
  ThermodynamicState state=setState_phX(p=p, h=h, X=X, phase=phase) "the state, evaluated once";
  HelmholtzCache cache=setHelmholtzCache(d=state.d, T=state.T, phase=state.phase) "evaluated once";
  Real d_der;
algorithm
  // in two-phase, T depends on p only (cp is infinite there)
  d_der := p_der*density_derp_h_cache(d=state.d, T=state.T, phase=state.phase, cache=cache) + h_der*density_derh_p_cache(d=state.d, T=state.T, phase=state.phase, cache=cache);
  state_der := ThermodynamicState(
    phase=0,
    T=p_der*jouleThomsonCoefficient_cache(d=state.d, T=state.T, phase=state.phase, cache=cache)
      + (if state.phase == 2 then 0 else h_der/specificHeatCapacityCp_cache(d=state.d, T=state.T, phase=state.phase, cache=cache)),
    p=p_der,
    d=d_der,
    u=h_der - p_der/state.d + p/state.d^2*d_der,
    h=h_der,
    s=(h_der - p_der/state.d)/state.T);
  annotation (Documentation(info="<html>
<p>
Derivative annotation of setState_phX. The output has the same record type as the state, with the time derivatives
of the Real components and phase=0. (The Modelica Language Specification 3.8-dev, section 12.7.1, describes a derivative
record with the Real components only; Dymola 2027x expects the same record type as the output of the original function.)
Not inlined: the state and the HelmholtzCache are evaluated once, and the partial derivatives are computed by the *_cache workers.
Dymola does not inline derivative functions inside the analytic Jacobians of nonlinear systems, where an inlined version evaluated the state several times per call.
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_phX_der;
