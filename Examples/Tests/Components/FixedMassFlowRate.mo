within HelmholtzMedia.Examples.Tests.Components;
model FixedMassFlowRate
  "Ideal pump that produces a constant mass flow rate from a large reservoir at fixed temperature and mass fraction"

  parameter Medium.MassFlowRate m_flow
    "Fixed mass flow rate from an infinite reservoir to the fluid port";

  parameter Boolean use_T_ambient=true "select T_ambient or h_ambient"
    annotation (Evaluate=true, Dialog(group=
          "Ambient temperature or ambient specific enthalpy"));
  parameter Medium.Temperature T_ambient=
      Modelica.SIunits.Conversions.from_degC(20) "Ambient temperature"
    annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                              enable=
          use_T_ambient));
  parameter Medium.SpecificEnthalpy h_ambient=
      1.e4 "Ambient specific enthalpy"
    annotation (Dialog(group="Ambient temperature or ambient specific enthalpy",
                                                              enable=not
          use_T_ambient));
  parameter Medium.MassFraction X_ambient[Medium.nX]
    "Ambient mass fractions m_i/m of reservoir";

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model"
     annotation (__Dymola_choicesAllMatching=true);

  Medium.BaseProperties medium "Medium in the source";
  FluidPort_b port(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
          rotation=0)));
equation
   if use_T_ambient then
     medium.T = T_ambient;
   else
     medium.h = h_ambient;
   end if;

   medium.Xi     = X_ambient[1:Medium.nXi];
   medium.p      = port.p;
   port.m_flow   = -m_flow;
   port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
   port.H_flow   = semiLinear(port.m_flow, port.h, medium.h);
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(
          extent={{20,60},{100,-60}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={192,192,192}),
        Rectangle(
          extent={{38,40},{100,-40}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255}),
        Ellipse(
          extent={{-100,80},{60,-80}},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          lineColor={0,0,255}),
        Polygon(
          points={{-60,70},{60,0},{-60,-68},{-60,70}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-54,32},{16,-30}},
          lineColor={255,0,0},
          textString="m"),
        Text(
          extent={{-142,142},{156,88}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-154,-88},{150,-132}},
          lineColor={0,0,0},
          textString="%m_flow"),
        Ellipse(
          extent={{-26,30},{-18,22}},
          lineColor={255,0,0},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>

</html>"));
end FixedMassFlowRate;
