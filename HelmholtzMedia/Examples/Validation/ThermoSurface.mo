within HelmholtzMedia.Examples.Validation;
model ThermoSurface "visualisation of thermodynamic surface"

  package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.ThermodynamicState state = Medium.setState_pT(p=p, T=T, phase=0);
  Medium.EoS.HelmholtzDerivs f = Medium.EoS.setHelmholtzDerivsThird(T=state.T, d=state.d, phase=1);

Modelica.Blocks.Sources.Ramp p_sub(
    height=pcrit - pmin,
    offset=pmin,
    duration=5,
    startTime=0)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));

Modelica.Blocks.Sources.Ramp p_super(
    height=pmax - pcrit,
    offset=0,
    duration=5,
    startTime=6)
    annotation (Placement(transformation(extent={{-40,20},{-20,40}})));

Modelica.Blocks.Sources.Ramp T_sub(
    height=Tcrit - Tmin,
    offset=Tmin,
    duration=5,
    startTime=0)
    annotation (Placement(transformation(extent={{-80,-40},{-60,-20}})));

Modelica.Blocks.Sources.Ramp T_super(
    height=Tmax - Tcrit,
    offset=0,
    duration=5,
    startTime=6)
    annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));

protected
  final constant Boolean appendToFile = false;
  final constant String fileName = "ThermoSurface.csv";
  final constant String Separator = ",";

  final constant Medium.Temperature Tmin=Medium.fluidLimits.TMIN;
  final constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  final constant Medium.Temperature Tmax=Medium.fluidLimits.TMAX;

  final constant Medium.AbsolutePressure pmin=1e-6;//Medium.fluidLimits.PMIN;
  final constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  final constant Medium.AbsolutePressure pmax=Medium.fluidLimits.PMAX;

  Medium.Temperature T(start=298.15);
  Medium.AbsolutePressure p_melt(start=101325) = Medium.Ancillary.meltingPressure_T(T);
  Medium.AbsolutePressure p(start=101325);

algorithm
  T := T_sub.y + T_super.y;
  p := min({p_melt, p_sub.y + p_super.y});

  if (time<=0) then
    if not appendToFile then
      // remove old file
      Modelica.Utilities.Files.remove(fileName);
    end if;
    // print headers
    Modelica.Utilities.Streams.print("d" +Separator
                                   + "T" +Separator
                                   + "p" +Separator
                                   + "s" +Separator
                                   + "u" +Separator
                                   + "h" +Separator
                                   + "g" +Separator
                                   + "dpdT" +Separator
                                   + "dpTd" +Separator
                                   + "dsdT" +Separator
                                   + "dudT" +Separator
                                   + "duTd" +Separator
                                   + "dhdT" +Separator
                                   + "dhTd" +Separator
                                   + "dgdT" +Separator
                                   + "dgTd" +Separator,
                                     fileName);
  end if;

  // print the actual values
  Modelica.Utilities.Streams.print(String(f.d) + Separator
                                 + String(f.T) + Separator
                                 + String(Medium.EoS.p(f)) + Separator
                                 + String(Medium.EoS.s(f)) + Separator
                                 + String(Medium.EoS.u(f)) + Separator
                                 + String(Medium.EoS.h(f))+Separator
                                 + String(Medium.EoS.g(f))+Separator
                                 + String(Medium.EoS.dpdT(f)) + Separator
                                 + String(Medium.EoS.dpTd(f)) + Separator
                                 + String(Medium.EoS.dsdT(f)) + Separator
                                 + String(Medium.EoS.dsTd(f)) + Separator
                                 + String(Medium.EoS.dudT(f)) + Separator
                                 + String(Medium.EoS.duTd(f)) + Separator
                                 + String(Medium.EoS.dhdT(f))+Separator
                                 + String(Medium.EoS.dhTd(f))+Separator
                                 + String(Medium.EoS.dgdT(f))+Separator
                                 + String(Medium.EoS.dgTd(f))+Separator,
                                   fileName);

annotation (experiment(StopTime=12, __Dymola_NumberOfIntervals=1000),
                                     __Dymola_experimentSetupOutput);
end ThermoSurface;
