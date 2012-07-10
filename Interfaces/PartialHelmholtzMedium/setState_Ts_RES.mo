within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_Ts_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: d=u) and output y (Residual)

  input Temperature T;
  input SpecificEntropy s;
  input FixedPhase phase=0;

protected
  Real R=Modelica.Constants.R/fluidConstants[1].molarMass
    "specific gas constant in J/kg.K";
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta=u/d_crit "reduced density";
  Real tau=T_crit/T "inverse reduced temperature";

  SpecificEntropy s_of_u;

algorithm
  s_of_u := R*(tau*(f_it(tau=tau, delta=delta) + f_rt(tau=tau, delta=delta)) - f_i(tau=tau, delta=delta) - f_r(tau=tau, delta=delta));
  // return the RESidual
  y := s - s_of_u;
end setState_Ts_RES;
