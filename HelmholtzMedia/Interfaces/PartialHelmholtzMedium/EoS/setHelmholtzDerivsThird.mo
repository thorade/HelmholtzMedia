within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivsThird

  input Density d;
  input Temperature T;
  input FixedPhase phase=1;
  output HelmholtzDerivs f;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";

algorithm
  f.d := d;
  f.T := T;
  f.delta := delta;
  f.tau := tau;

  f.i    := f_i(tau=tau, delta=delta);
  f.it   := f_it(tau=tau, delta=delta);
  f.itt  := f_itt(tau=tau, delta=delta);
  f.ittt := f_ittt(tau=tau, delta=delta);

  f.r    := f_r(tau=tau, delta=delta);
  f.rt   := f_rt(tau=tau, delta=delta);
  f.rtt  := f_rtt(tau=tau, delta=delta);
  f.rttt := f_rttt(tau=tau, delta=delta);
  f.rtd  := f_rtd(tau=tau, delta=delta);
  f.rttd := f_rttd(tau=tau, delta=delta);
  f.rtdd := f_rtdd(tau=tau, delta=delta);
  f.rd   := f_rd(tau=tau, delta=delta);
  f.rdd  := f_rdd(tau=tau, delta=delta);
  f.rddd := f_rddd(tau=tau, delta=delta);

  assert(phase==1, "This function will return valid values for single phase input only!");

end setHelmholtzDerivsThird;
