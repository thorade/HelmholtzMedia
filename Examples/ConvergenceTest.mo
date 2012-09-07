within HelmholtzMedia.Examples;
package ConvergenceTest

  model setSat
    package medium = HelmholtzFluids.Butane;
    medium.SaturationProperties sat;
    medium.SaturationProperties sat_p;
    medium.SaturationProperties sat_dl;
    medium.SaturationProperties sat_dv;

    Modelica.Blocks.Sources.Ramp T_ramp(
      duration=8,
      startTime=1,
      height=Tcrit - Tmin,
      offset=Tmin)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

  protected
    constant medium.Temperature Tmin=medium.fluidLimits.TMIN;
    constant medium.Temperature Tcrit=medium.fluidConstants[1].criticalTemperature;
    constant medium.Temperature Tmax=medium.fluidLimits.TMAX;

  equation
    // forward
    sat = medium.setSat_T(T=T_ramp.y);

    // backward
    sat_p = medium.setSat_p(p=sat.psat);
    sat_dl = medium.setSat_d(d=sat.liq.d);
    sat_dv = medium.setSat_d(d=sat.vap.d);

    annotation (experiment(
        StopTime=10,
        NumberOfIntervals=10000,
        Tolerance=1e-005));
  end setSat;

  model SinglePhase_setState
    package Medium = HelmholtzFluids.Butane;
    Medium.AbsolutePressure p;
    Medium.Temperature T;
    Medium.ThermodynamicState state;
    Medium.ThermodynamicState state_dT;
    Medium.ThermodynamicState state_pd;
    Medium.ThermodynamicState state_ph;
    Medium.ThermodynamicState state_ps;
    Medium.ThermodynamicState state_Ts;

  Modelica.Blocks.Sources.Ramp T_sub(
      duration=4,
      startTime=0.1,
      height=Tcrit - Tmin,
      offset=Tmin)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

  Modelica.Blocks.Sources.Ramp T_super(
      height=Tmax - Tcrit,
      duration=5,
      startTime=6,
      offset=0)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

  Modelica.Blocks.Sources.Sine p_sine(
      amplitude=(pmax - pmin)/2,
      offset=(pmax - pmin)/2,
      freqHz=100,
      startTime=0)
      annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

  protected
    constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
    constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
    constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
    constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
    constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
    constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

  equation
    T = T_sub.y + T_super.y;
    p = p_sine.y;

    // set state to valid single-phase values
    state=Medium.setState_pT(p=p, T=T, phase=0);

    // call other setState functions
    state_dT=Medium.setState_dT(d=state.d, T=state.T, phase=0);
    state_pd=Medium.setState_pd(p=state.p, d=state.d, phase=0);
    state_ph=Medium.setState_ph(p=state.p, h=state.h, phase=0);
    state_ps=Medium.setState_ps(p=state.p, s=state.s, phase=0);
    state_Ts=Medium.setState_Ts(T=state.T, s=state.s, phase=0);

    annotation (experiment(StopTime=12, NumberOfIntervals=10000));
  end SinglePhase_setState;

  model TwoPhase_setState
    package Medium = HelmholtzFluids.Butane;
    Medium.Density d;
    Medium.Temperature T;
    Medium.ThermodynamicState state;
    Medium.ThermodynamicState state_dT;
    Medium.ThermodynamicState state_pd;
    Medium.ThermodynamicState state_ph;
    Medium.ThermodynamicState state_ps;
    Medium.ThermodynamicState state_Ts;

  Modelica.Blocks.Sources.Ramp T_sub(
      height=Tcrit - Tmin,
      offset=Tmin,
      duration=9,
      startTime=1)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

  protected
    constant Medium.Density dmin=Medium.fluidLimits.DMIN;
    constant Medium.Density dcrit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
    constant Medium.Density dmax=Medium.fluidLimits.DMAX;
    constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
    constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
    constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;

  equation
    d = dcrit;
    T = T_sub.y;

    // set state to valid two-phase values
    state=Medium.setState_Tx(T=T, x=0.5);

    // call setState functions
    state_dT=Medium.setState_dT(d=state.d, T=state.T, phase=0);
    state_pd=Medium.setState_pd(p=state.p, d=state.d, phase=0);
    state_ph=Medium.setState_ph(p=state.p, h=state.h, phase=0);
    state_ps=Medium.setState_ps(p=state.p, s=state.s, phase=0);
    state_Ts=Medium.setState_Ts(T=state.T, s=state.s, phase=0);

    annotation (experiment(StopTime=12, NumberOfIntervals=10000));
  end TwoPhase_setState;

  model AncillaryFunctions
    package medium = HelmholtzFluids.Butane;
    medium.Temperature Tsat;
    medium.AbsolutePressure psat;
    medium.Density dliq;
    medium.Density dvap;

    medium.Temperature T_p;
    medium.Temperature T_dl;
    medium.Temperature T_dv;

    Modelica.Blocks.Sources.Ramp T_ramp(
      duration=5,
      height=T_crit - T_trip,
      offset=T_trip,
      startTime=1)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

  protected
    constant medium.Temperature T_trip=medium.fluidConstants[1].triplePointTemperature;
    constant medium.Temperature T_crit=medium.fluidConstants[1].criticalTemperature;

  algorithm
    // forward
    Tsat := T_ramp.y;
    psat := medium.Ancillary.saturationPressure_T(T=Tsat);
    dliq := medium.Ancillary.bubbleDensity_T(T=Tsat);
    dvap := medium.Ancillary.dewDensity_T(T=Tsat);

    // inverse
    T_p  := medium.Ancillary.saturationTemperature_p(p=psat);
    T_dl := medium.Ancillary.saturationTemperature_d(d=dliq);
    T_dv := medium.Ancillary.saturationTemperature_d(d=dvap);

  annotation (experiment(
        StopTime=11,
        NumberOfIntervals=1000,
        Tolerance=1e-005), __Dymola_experimentSetupOutput);
  end AncillaryFunctions;

  model SinglePhase_Transport
    package Medium = HelmholtzFluids.Butane;
    Medium.AbsolutePressure p;
    Medium.Temperature T;
    Medium.ThermodynamicState state;
    Medium.DynamicViscosity eta;
    Medium.ThermalConductivity lambda;

  Modelica.Blocks.Sources.Ramp T_sub(
      height=Tcrit - Tmin,
      duration=4,
      offset=Tmin,
      startTime=0.1)
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.Blocks.Sources.Sine p_sine(
      phase=0,
      amplitude=(pmax - pmin)/2,
      offset=(pmax - pmin)/2,
      startTime=0,
      freqHz=100)
      annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

  Modelica.Blocks.Sources.Ramp T_super(
      height=Tmax - Tcrit,
      duration=5,
      startTime=6,
      offset=0)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

  protected
    constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
    constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
    constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;
    constant Medium.AbsolutePressure pmin=Medium.fluidLimits.PMIN;
    constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
    constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

  equation
    T = T_sub.y + T_super.y;
    p = p_sine.y;

    state=Medium.setState_pT(p=p, T=T, phase=0);
    eta=Medium.dynamicViscosity(state);
    lambda=Medium.thermalConductivity(state);

    annotation (experiment(StopTime=12, NumberOfIntervals=10000));
  end SinglePhase_Transport;
end ConvergenceTest;
