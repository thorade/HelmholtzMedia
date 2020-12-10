within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivsZero

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
  f.r   := f_r(tau=f.tau, delta=f.delta);

  assert(phase==1, "This function will return valid values for single phase input only!");

end setHelmholtzDerivsZero;
