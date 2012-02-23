within HelmholtzFluids.Examples.Tests.Components;
model ShortPipe "Simple pressure loss in pipe"
   replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model"
     annotation (__Dymola_choicesAllMatching=true);

  parameter Medium.AbsolutePressure dp_nominal(min=1.e-10)
    "Nominal pressure drop";
  parameter Medium.MassFlowRate m_flow_nominal(min=1.e-10)
    "Nominal mass flow rate at nominal pressure drop";

  FluidPort_a port_a(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-120,-10},{-100,10}},
          rotation=0)));
  FluidPort_b port_b(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{120,-10},{100,10}},
          rotation=0)));
  // Medium.BaseProperties medium_a(p=port_a.p, h=port_a.h, Xi=port_a.Xi)
  //   "Medium properties in port_a";
  // Medium.BaseProperties medium_b(p=port_b.p, h=port_b.h, Xi=port_b.Xi)
  //   "Medium properties in port_b";
  Medium.MassFlowRate m_flow
    "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
  Modelica.SIunits.Pressure dp "Pressure drop from port_a to port_b";
equation
  /* Handle reverse and zero flow */
  port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,   port_b.h);
  port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);

  /* Energy, mass and substance mass balance */
  port_a.H_flow + port_b.H_flow = 0;
  port_a.m_flow + port_b.m_flow = 0;
  port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);

  // Design direction of mass flow rate
  m_flow = port_a.m_flow;

  // Pressure drop
  dp = port_a.p - port_b.p;
  m_flow = (m_flow_nominal/dp_nominal)*dp;
                                                                                   annotation (Icon(
        coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Rectangle(
          extent={{-100,60},{100,-60}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={192,192,192}),
        Rectangle(
          extent={{-100,34},{100,-36}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255}),
        Text(
          extent={{-150,140},{150,80}},
          lineColor={0,0,0},
          textString="%name"),
        Text(
          extent={{-136,-62},{122,-108}},
          lineColor={0,0,0},
          textString="k=%m_flow_nominal/%dp_nominal")}),
                                           Documentation(info="<html>
<p>
Model <b>ShortPipe</b> defines a simple pipe model
with pressure loss due to friction. It is assumed that
no mass or energy is stored in the pipe.
</p>
</html>"));
end ShortPipe;
