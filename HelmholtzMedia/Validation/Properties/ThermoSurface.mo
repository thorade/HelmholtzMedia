within HelmholtzMedia.Validation.Properties;
model ThermoSurface "visualisation of thermodynamic surface"

  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.ThermodynamicState state = Medium.setState_pT(p=p, T=T, phase=0);
  Medium.Internal.EoS.HelmholtzDerivs f = Medium.Internal.EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=1, order=3);

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
  Medium.AbsolutePressure p_melt(start=101325) = Medium.Internal.Ancillary.meltingPressure_T(T);
  Medium.AbsolutePressure p(start=101325);

equation
  T = T_sub.y + T_super.y;
  p = min(p_melt, p_sub.y + p_super.y);

algorithm
  // only file output here; T and p are computed in equations above to avoid an algebraic loop
  when initial() then
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
  end when;

  // print the actual values
  Modelica.Utilities.Streams.print(String(f.d) + Separator
                                 + String(f.T) + Separator
                                 + String(Medium.Internal.EoS.p(f)) + Separator
                                 + String(Medium.Internal.EoS.s(f)) + Separator
                                 + String(Medium.Internal.EoS.u(f)) + Separator
                                 + String(Medium.Internal.EoS.h(f))+Separator
                                 + String(Medium.Internal.EoS.g(f))+Separator
                                 + String(Medium.Internal.EoS.dpdT(f)) + Separator
                                 + String(Medium.Internal.EoS.dpTd(f)) + Separator
                                 + String(Medium.Internal.EoS.dsdT(f)) + Separator
                                 + String(Medium.Internal.EoS.dsTd(f)) + Separator
                                 + String(Medium.Internal.EoS.dudT(f)) + Separator
                                 + String(Medium.Internal.EoS.duTd(f)) + Separator
                                 + String(Medium.Internal.EoS.dhdT(f))+Separator
                                 + String(Medium.Internal.EoS.dhTd(f))+Separator
                                 + String(Medium.Internal.EoS.dgdT(f))+Separator
                                 + String(Medium.Internal.EoS.dgTd(f))+Separator,
                                   fileName);

  annotation (experiment(StopTime=12));
end ThermoSurface;
