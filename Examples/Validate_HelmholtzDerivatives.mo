within HelmholtzMedia.Examples;
model Validate_HelmholtzDerivatives
  package medium = HelmholtzFluids.Butane;
  Real tau;
  Real delta;

  Real ai;
  Real ai_tau;
  Real ai_tau_tau;

  Real ar;
  Real ar_delta;
  Real ar_delta_delta;
  Real ar_delta_tau;
  Real ar_tau;
  Real ar_tau_tau;

  Modelica.Blocks.Sources.Ramp tau_ramp(
    height=2,
    offset=Modelica.Constants.eps,
    startTime=1,
    duration=8)
    annotation (Placement(transformation(extent={{-60,22},{-40,42}})));
  Modelica.Blocks.Sources.Ramp delta_ramp(
    height=2,
    offset=Modelica.Constants.eps,
    duration=8,
    startTime=1)
              annotation (Placement(transformation(extent={{-58,-20},{-38,0}})));

algorithm
    tau := tau_ramp.y;
    delta := delta_ramp.y;

    ai :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_i(delta=delta, tau=tau);
    ai_tau :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_it(delta=delta, tau=tau);
    ai_tau_tau :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_itt(delta=delta, tau=tau);

    ar :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_r(delta=delta, tau=tau);
    ar_delta :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_rd(delta=delta, tau=tau);
    ar_delta_delta :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_rdd(delta=delta, tau=tau);
    ar_delta_tau :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_rtd(delta=delta, tau=tau);
    ar_tau :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_rt(delta=delta, tau=tau);
    ar_tau_tau :=HelmholtzFluids.Interfaces.PartialHelmholtzMedium.a_rtt(delta=delta, tau=tau);

  annotation (experiment(StopTime=10), __Dymola_experimentSetupOutput);
end Validate_HelmholtzDerivatives;
