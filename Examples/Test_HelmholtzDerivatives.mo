within HelmholtzFluids.Examples;
model Test_HelmholtzDerivatives
  package medium = HelmholtzFluids.R134a;
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
    duration=6,
    height=2,
    offset=Modelica.Constants.eps)
    annotation (Placement(transformation(extent={{-60,22},{-40,42}})));
  Modelica.Blocks.Sources.Ramp delta_ramp(
    duration=7,
    height=2,
    offset=Modelica.Constants.eps)
              annotation (Placement(transformation(extent={{-58,-20},{-38,0}})));

algorithm
    tau := tau_ramp.y;
    delta := delta_ramp.y;

    ai :=medium.ai(delta=delta, tau=tau);
    ai_tau :=medium.ai_tau(delta=delta, tau=tau);
    ai_tau_tau :=medium.ai_tau_tau(delta=delta, tau=tau);

    ar :=medium.ar(delta=delta, tau=tau);
    ar_delta :=medium.ar_delta(delta=delta, tau=tau);
    ar_delta_delta :=medium.ar_delta_delta(delta=delta, tau=tau);
    ar_delta_tau :=medium.ar_delta_tau(delta=delta, tau=tau);
    ar_tau :=medium.ar_tau(delta=delta, tau=tau);
    ar_tau_tau :=medium.ar_tau_tau(delta=delta, tau=tau);

  annotation (experiment(StopTime=10), __Dymola_experimentSetupOutput);
end Test_HelmholtzDerivatives;
