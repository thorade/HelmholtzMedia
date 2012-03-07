within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setHelmholtzDerivs

  input Temperature T;
  input Density d;
  input FixedPhase phase(min=0, max=2);
  output HelmholtzDerivs f;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";

algorithm
  f.T := T;
  f.d := d;
  f.R := R;
  f.tau := tau;
  f.delta := delta;

  if (phase==1) then
    f.i   := f_i(tau=tau, delta=delta);
    f.it  := f_it(tau=tau, delta=delta);
    f.itt := f_itt(tau=tau, delta=delta);

    f.r   := f_r(tau=tau, delta=delta);
    f.rt  := f_rt(tau=tau, delta=delta);
    f.rtt := f_rtt(tau=tau, delta=delta);
    f.rtd := f_rtd(tau=tau, delta=delta);
    f.rd  := f_rd(tau=tau, delta=delta);
    f.rdd := f_rdd(tau=tau, delta=delta);
  // else: do nothing
  end if;
end setHelmholtzDerivs;
