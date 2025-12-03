within HelmholtzMedia.Examples.MediaTestModels;
partial model PartialTestModel "Basic test model to test a medium"
  import SI = Modelica.SIunits;

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model" annotation (choicesAllMatching=true);
  parameter SI.AbsolutePressure p_start=Medium.p_default "Initial value of pressure";
  parameter SI.Temperature T_start=Medium.T_default "Initial value of temperature";
  parameter SI.SpecificEnthalpy h_start=Medium.h_default "Initial value of specific enthalpy";
  parameter Real X_start[Medium.nX]=Medium.X_default "Initial value of mass fractions";

  Modelica.Media.Examples.Tests.Components.PortVolume volume(
    redeclare package Medium = Medium,
    p_start=p_start,
    T_start=T_start,
    h_start=h_start,
    X_start=X_start,
    V=0.1) annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Modelica.Media.Examples.Tests.Components.FixedMassFlowRate fixedMassFlowRate(
    redeclare package Medium = Medium,
    T_ambient=1.2*T_start,
    h_ambient=1.2*h_start,
    m_flow=1,
    X_ambient=0.5*X_start) annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
  Modelica.Media.Examples.Tests.Components.FixedAmbient ambient(
    redeclare package Medium = Medium,
    T_ambient=T_start,
    h_ambient=h_start,
    X_ambient=X_start,
    p_ambient=p_start) annotation (Placement(transformation(extent={{60,0},{40,20}})));
  Modelica.Media.Examples.Tests.Components.ShortPipe shortPipe(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=0.1e5) annotation (Placement(transformation(extent={{0,0},{20,20}})));
equation
  connect(fixedMassFlowRate.port, volume.port) annotation (Line(points={{-59,10},{-30,10}}, color={0,127,255}));
  connect(volume.port, shortPipe.port_a) annotation (Line(points={{-30,10},{-1,10}}, color={0,127,255}));
  connect(shortPipe.port_b, ambient.port) annotation (Line(points={{21,10},{39,10}}, color={0,127,255}));
  annotation (Documentation(info="<html>

</html>"), obsolete="Obsolete model - use Modelica.Media.Examples.Utilities.PartialTestModel instead");
end PartialTestModel;
