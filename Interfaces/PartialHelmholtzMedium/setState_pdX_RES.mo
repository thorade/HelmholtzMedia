within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_pdX_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

  input AbsolutePressure p;
  input Density d;
  input FixedPhase phase=0;

protected
  Real R=Modelica.Constants.R/fluidConstants[1].molarMass
    "specific gas constant in J/kg.K";
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta=d/d_crit "reduced density";
  Real tau=T_crit/u "inverse reduced temperature";

  AbsolutePressure p_of_u;

algorithm
  p_of_u := ((1 + delta*f_rd(delta=delta, tau=tau))*d*R*u);
  y := p - p_of_u;
end setState_pdX_RES;
