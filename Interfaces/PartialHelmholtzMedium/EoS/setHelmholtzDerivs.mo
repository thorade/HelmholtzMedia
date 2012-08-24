within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function setHelmholtzDerivs

  input Density d;
  input Temperature T;
  input FixedPhase phase=1;
  output HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.HelmholtzDerivs
                         f;

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
    f.i   := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_i(
                 tau=tau, delta=delta);
    f.it  := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_it(
                  tau=tau, delta=delta);
    f.itt := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_itt(
                   tau=tau, delta=delta);

    f.r   := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_r(
                 tau=tau, delta=delta);
    f.rt  := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_rt(
                  tau=tau, delta=delta);
    f.rtt := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_rtt(
                   tau=tau, delta=delta);
    f.rtd := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_rtd(
                   tau=tau, delta=delta);
    f.rd  := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_rd(
                  tau=tau, delta=delta);
    f.rdd := HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.f_rdd(
                   tau=tau, delta=delta);
  // else: do nothing
  end if;
end setHelmholtzDerivs;
