within HelmholtzFluids.Examples.Tests.Components;
partial model PartialTestModel2 "slightly larger test model to test a medium"
  import SI = Modelica.SIunits;

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model" annotation (__Dymola_choicesAllMatching=true);
  parameter SI.AbsolutePressure p_start = 1.0e5 "Initial value of pressure";
  parameter SI.Temperature T_start = 300 "Initial value of temperature";
  parameter SI.SpecificEnthalpy h_start = 1
    "Initial value of specific enthalpy";
  parameter Real X_start[Medium.nX] = Medium.reference_X
    "Initial value of mass fractions";
  PortVolume volume(redeclare package Medium = Medium,
                    p_start=p_start,
                    T_start=T_start,
                    h_start=h_start,
                    X_start = X_start,
                    V=0.1)
           annotation (Placement(transformation(extent={{-60,0},{-40,20}},
          rotation=0)));
  FixedMassFlowRate fixedMassFlowRate(redeclare package Medium = Medium,
    T_ambient=1.2*T_start,
    h_ambient=1.2*h_start,
    m_flow=1,
    X_ambient=0.5*X_start)
                          annotation (Placement(transformation(extent={{
            -100,0},{-80,20}}, rotation=0)));
  FixedAmbient ambient(
    redeclare package Medium = Medium,
    T_ambient=T_start,
    h_ambient=h_start,
    X_ambient=X_start,
    p_ambient=p_start) annotation (Placement(transformation(extent={{92,0},
            {72,20}}, rotation=0)));
  ShortPipe shortPipe(redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=0.1e5)
    annotation (Placement(transformation(extent={{-30,0},{-10,20}},
          rotation=0)));
  PortVolume volume1(
                    redeclare package Medium = Medium,
                    p_start=p_start,
                    T_start=T_start,
                    h_start=h_start,
                    X_start = X_start,
                    V=0.1)
           annotation (Placement(transformation(extent={{0,0},{20,20}},
          rotation=0)));
  ShortPipe shortPipe1(
                      redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=0.1e5)
    annotation (Placement(transformation(extent={{36,0},{56,20}},
          rotation=0)));
equation
  connect(fixedMassFlowRate.port, volume.port) annotation (Line(points={{
          -79,10},{-50,10}}, color={0,127,255}));
  connect(volume.port, shortPipe.port_a)
    annotation (Line(points={{-50,10},{-31,10}}, color={0,127,255}));
  connect(volume1.port, shortPipe1.port_a)
    annotation (Line(points={{10,10},{35,10}}, color={0,127,255}));
  connect(shortPipe.port_b, volume1.port)
    annotation (Line(points={{-9,10},{10,10}}, color={0,127,255}));
  connect(shortPipe1.port_b, ambient.port)
    annotation (Line(points={{57,10},{71,10}}, color={0,127,255}));
  annotation (         Documentation(info="<html>

</html>"));
end PartialTestModel2;
