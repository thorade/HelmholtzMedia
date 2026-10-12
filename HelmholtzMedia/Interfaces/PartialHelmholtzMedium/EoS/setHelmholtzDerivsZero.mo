within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivsZero

  input Density d;
  input Temperature T;
  input FixedPhase phase=1;
  output HelmholtzDerivs f;

algorithm
  f.d := d;
  f.T := T;
  f.delta := d/f.d_reducing;
  f.tau := f.T_reducing/T;

  (f.i, f.it, f.itt, f.ittt) := f_i_derivs(delta=f.delta, tau=f.tau, order=0);
  (f.r, f.rd, f.rdd, f.rt, f.rtt, f.rtd, f.rttt, f.rttd, f.rtdd, f.rddd) := f_r_derivs(delta=f.delta, tau=f.tau, order=0);

  assert(phase==1, "This function will return valid values for single phase input only!");

end setHelmholtzDerivsZero;
