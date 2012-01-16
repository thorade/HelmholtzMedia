within HelmholtzFluids.PartialHelmholtzFluid;
record EosLimits "validity limits for fluid model"
  extends Modelica.Icons.Record;
  constant Modelica.SIunits.Temperature TMIN "minimum temperature";
  constant Modelica.SIunits.Temperature TNOM "nominal temperature";
  constant Modelica.SIunits.Temperature TMAX "maximum temperature";
  constant Modelica.SIunits.Density DMIN "minimum density";
  constant Modelica.SIunits.Density DNOM "nominal density";
  constant Modelica.SIunits.Density DMAX "maximum density";
  constant Modelica.SIunits.AbsolutePressure PMIN "minimum pressure";
  constant Modelica.SIunits.AbsolutePressure PNOM "nominal pressure";
  constant Modelica.SIunits.AbsolutePressure PMAX "maximum pressure";
  constant Modelica.SIunits.SpecificEnthalpy HMIN "minimum enthalpy";
  constant Modelica.SIunits.SpecificEnthalpy HMAX "maximum enthalpy";
  constant Modelica.SIunits.SpecificEntropy SMIN "minimum entropy";
  constant Modelica.SIunits.SpecificEntropy SMAX "maximum entropy";
end EosLimits;
