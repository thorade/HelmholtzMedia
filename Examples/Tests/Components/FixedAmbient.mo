within HelmholtzFluids.Examples.Tests.Components;
model FixedAmbient "Ambient pressure, temperature and mass fraction source"
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model"
     annotation (__Dymola_choicesAllMatching=true);

  parameter Boolean use_p_ambient=true "select p_ambient or d_ambient"
    annotation (Evaluate=true, Dialog(group=
          "Ambient pressure or ambient density"));
  parameter Medium.AbsolutePressure p_ambient= 101325 "Ambient pressure"          annotation (
     Dialog(group="Ambient pressure or ambient density", enable=use_p_ambient));
  parameter Medium.Density d_ambient=1 "Ambient density"
                       annotation (Dialog(group=
          "Ambient pressure or ambient density", enable=not use_p_ambient));
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
    "Ambient mass fractions m_i/m"                                                   annotation (Dialog(group=
          "Only for multi-substance flow", enable=Medium.nX > 0));

  Medium.BaseProperties medium "Medium in the source";
  FluidPort_b port(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
          rotation=0)));

equation
  if use_p_ambient or Medium.singleState then
    medium.p = p_ambient;
  else
    medium.d = d_ambient;
  end if;

  if use_T_ambient then
    medium.T = T_ambient;
  else
    medium.h = h_ambient;
  end if;

  medium.Xi = X_ambient[1:Medium.nXi];

  port.p = medium.p;
  port.H_flow   = semiLinear(port.m_flow, port.h, medium.h);
  port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Ellipse(
          extent={{-100,80},{100,-80}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Sphere,
          fillColor={0,127,255}), Text(
          extent={{-136,144},{132,82}},
          textString="%name",
          lineColor={0,0,255})}),
    Documentation(info="<html>
<p>
Model <b>FixedAmbient_pt</b> defines constant values for ambient conditions:
</p>
<ul>
<li> Ambient pressure.</li>
<li> Ambient temperature.</li>
<li> Ambient mass fractions (only for multi-substance flow).</li>
</ul>
<p>
Note, that ambient temperature
and mass fractions have only an effect if the mass flow
is from the ambient into the port. If mass is flowing from
the port into the ambient, the ambient definitions,
with exception of ambient pressure, do not have an effect.
</p>
</html>"));
end FixedAmbient;
