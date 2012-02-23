within HelmholtzFluids.Examples.Tests.Components;
model PortVolume
  "Fixed volume associated with a port by the finite volume method"
  import SI = Modelica.SIunits;

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model"
     annotation (__Dymola_choicesAllMatching=true);

  parameter SI.Volume V=1e-6 "Fixed size of junction volume";

  parameter Boolean use_p_start=true "select p_start or d_start"
    annotation (Evaluate=true, Dialog(group="Initial pressure or initial density"));
  parameter Medium.AbsolutePressure p_start = 101325 "Initial pressure"
    annotation (Dialog(group="Initial pressure or initial density", enable=use_p_start));
  parameter Medium.Density d_start=1 "Initial density"
    annotation (Dialog(group="Initial pressure or initial density", enable=not use_p_start));
  parameter Boolean use_T_start=true "select T_start or h_start"
    annotation (Evaluate=true, Dialog(group="Initial temperature or initial specific enthalpy"));
  parameter Medium.Temperature T_start = Modelica.SIunits.Conversions.from_degC(20)
    "Initial temperature"
    annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=use_T_start));
  parameter Medium.SpecificEnthalpy h_start = 1.e4 "Initial specific enthalpy"
    annotation (Dialog(group="Initial temperature or initial specific enthalpy", enable=not use_T_start));
  parameter Medium.MassFraction X_start[Medium.nX]
    "Initial mass fractions m_i/m"
    annotation (Dialog(group="Only for multi-substance flow", enable=Medium.nX > 0));

  FluidPort_a port(redeclare package Medium = Medium) annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, rotation=0)));
  Medium.BaseProperties medium(preferredMediumStates=true);
  SI.Energy U "Internal energy of port volume";
  SI.Mass m "Mass of junction volume";
  SI.Mass mXi[Medium.nXi] "Independent substance masses of junction volume";

initial equation
  if not Medium.singleState then
    if use_p_start then
       medium.p = p_start;
    else
       medium.d = d_start;
    end if;
  end if;

  if use_T_start then
     medium.T = T_start;
  else
     medium.h = h_start;
  end if;

  medium.Xi = X_start[1:Medium.nXi];
equation
  // Connect port to medium variables
     medium.p = port.p;
     medium.h = port.h;
     medium.Xi = port.Xi;

  // Total quantities
     m    = V*medium.d;
     mXi = m*medium.Xi;
     U    = m*medium.u;

  // Mass and energy balance
     der(m)    = port.m_flow;
     der(mXi) = port.mXi_flow;
     der(U)    = port.H_flow;
  annotation (
   Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Sphere,
          fillColor={170,213,255}),
        Text(
          extent={{-144,178},{146,116}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-130,-108},{144,-150}},
          lineColor={0,0,0},
          textString="V=%V")}),
                         Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport.
</p>
</html>"),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}),
            graphics));
end PortVolume;
