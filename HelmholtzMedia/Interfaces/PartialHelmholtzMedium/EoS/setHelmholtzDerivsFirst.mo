within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivsFirst

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
  f.r   := f_r(tau=f.tau, delta=f.delta);
  f.rt  := f_rt(tau=f.tau, delta=f.delta);
  f.rd  := f_rd(tau=f.tau, delta=f.delta);

  assert(phase==1, "This function will return valid values for single phase input only!");

end setHelmholtzDerivsFirst;
