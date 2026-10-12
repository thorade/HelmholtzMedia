within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivs
  "dimensionless Helmholtz energy and its derivatives up to the given order (single-phase input)"
  input Density d;
  input Temperature T;
  input Integer order(min=0, max=3)=2 "highest order of the derivatives to compute";
  input Boolean includeIdeal=true
    "false: skip the ideal part (fields i, it, itt, ittt are set to zero); sufficient for p and its derivatives";
  input FixedPhase phase=1;
  output HelmholtzDerivs f;

algorithm
  f.d := d;
  f.T := T;
  f.delta := d/f.d_reducing;
  f.tau := f.T_reducing/T;

  if includeIdeal then
    (f.i, f.it, f.itt, f.ittt) := f_i_derivs(delta=f.delta, tau=f.tau, order=order);
  else
    f.i := 0;
    f.it := 0;
    f.itt := 0;
    f.ittt := 0;
  end if;
  (f.r, f.rd, f.rdd, f.rt, f.rtt, f.rtd, f.rttt, f.rttd, f.rtdd, f.rddd) := f_r_derivs(delta=f.delta, tau=f.tau, order=order);

  assert(phase==1, "This function will return valid values for single phase input only!");

  annotation (Documentation(info="<html>
<p>
Fills the record HelmholtzDerivs for given density and temperature with the ideal and residual parts of the
dimensionless Helmholtz energy and their derivatives up to <code>order</code> (0 to 3). Derivatives of higher
order than requested are set to zero.
</p>
<p>
With <code>includeIdeal=false</code> the ideal part is not computed and its fields are set to zero.
This is sufficient for the pressure and its derivatives (EoS.p, EoS.dpdT, EoS.dpTd), e.g. inside the
iterations of setState_pTX and setState_pd, but not for h, s, u, g, cp, cv or w.
</p>
</html>"));
end setHelmholtzDerivs;
