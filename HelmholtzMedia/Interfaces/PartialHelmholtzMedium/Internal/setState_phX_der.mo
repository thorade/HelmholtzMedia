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
algorithm
  state_der := ThermodynamicState(
    phase=0,
    T=temperature_ph_state_der(p=p, h=h, state=setState_phX(p=p, h=h, X=X, phase=phase), p_der=p_der, h_der=h_der),
    p=p_der,
    d=density_ph_state_der(p=p, h=h, state=setState_phX(p=p, h=h, X=X, phase=phase), p_der=p_der, h_der=h_der),
    u=h_der - p_der/density(state=setState_phX(p=p, h=h, X=X, phase=phase))
      + p/density(state=setState_phX(p=p, h=h, X=X, phase=phase))^2*density_ph_state_der(p=p, h=h, state=setState_phX(p=p, h=h, X=X, phase=phase), p_der=p_der, h_der=h_der),
    h=h_der,
    s=specificEntropy_ph_state_der(p=p, h=h, state=setState_phX(p=p, h=h, X=X, phase=phase), p_der=p_der, h_der=h_der));
  annotation (Inline=true, Documentation(info="<html>
<p>
Derivative annotation of setState_phX. The output has the same record type as the state, with the time derivatives
of the Real components and phase=0. (The Modelica Language Specification 3.8-dev, section 12.7.1, describes a derivative
record with the Real components only; Dymola 2027x expects the same record type as the output of the original function.)
Written as one statement and inlined, so that the setState_phX calls are shared with the original call (same arguments).
u = h - p/d, so du/dt = dh/dt - (dp/dt)/d + p/d^2*(dd/dt).
</p>
</html>"));
end setState_phX_der;
