within HelmholtzMedia.Examples.Tests.Components;
connector FluidPort
  "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"

  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model" annotation (__Dymola_choicesAllMatching=true);

  Medium.AbsolutePressure p "Pressure in the connection point";
  flow Medium.MassFlowRate m_flow
    "Mass flow rate from the connection point into the component";

  Medium.SpecificEnthalpy h "Specific mixture enthalpy in the connection point";
  flow Medium.EnthalpyFlowRate H_flow
    "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)";

  Medium.MassFraction Xi[Medium.nXi]
    "Independent mixture mass fractions m_i/m in the connection point";
  flow Medium.MassFlowRate mXi_flow[Medium.nXi]
    "Mass flow rates of the independent substances from the connection point into the component (if m_flow > 0, mX_flow = m_flow*X)";

  Medium.ExtraProperty C[Medium.nC] "properties c_i/m in the connection point";
  flow Medium.ExtraPropertyFlowRate mC_flow[Medium.nC]
    "Flow rates of auxiliary properties from the connection point into the component (if m_flow > 0, mC_flow = m_flow*C)";

  annotation (Documentation(info="<html>

</html>"));
end FluidPort;
