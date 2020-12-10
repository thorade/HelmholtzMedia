within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivsSecond

  input Density d;
  input Temperature T;
  input FixedPhase phase=1;
  output HelmholtzDerivs f;

algorithm
  f.d := d;
  f.T := T;
  f.delta := d/f.d_crit;
  f.tau := f.T_crit/T;

  f.i   := f_i(tau=f.tau, delta=f.delta);
  f.it  := f_it(tau=f.tau, delta=f.delta);
  f.itt := f_itt(tau=f.tau, delta=f.delta);

  f.r   := f_r(tau=f.tau, delta=f.delta);
  f.rt  := f_rt(tau=f.tau, delta=f.delta);
  f.rtt := f_rtt(tau=f.tau, delta=f.delta);
  f.rtd := f_rtd(tau=f.tau, delta=f.delta);
  f.rd  := f_rd(tau=f.tau, delta=f.delta);
  f.rdd := f_rdd(tau=f.tau, delta=f.delta);

  assert(phase==1, "This function will return valid values for single phase input only!");

end setHelmholtzDerivsSecond;
