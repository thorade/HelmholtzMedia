within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_pTX_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: d=u) and output y (Residual)

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.AbsolutePressure p;

protected
  Real R=Modelica.Constants.R/fluidConstants[1].molarMass
    "specific gas constant in J/kg.K";
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta=u/d_crit "reduced density";
  Real tau=T_crit/T "inverse reduced temperature";

  AbsolutePressure p_of_u;

algorithm
  p_of_u := ((1 + delta*f_rd(delta=delta, tau=tau))*u*R*T);
  y := p - p_of_u;
end setState_pTX_RES;
