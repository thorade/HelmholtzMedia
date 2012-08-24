within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_ThX_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: d=u) and output y (Residual)

  input Temperature T;
  input SpecificEnthalpy h;
  input FixedPhase phase=0;

protected
  Real R=Modelica.Constants.R/fluidConstants[1].molarMass
    "specific gas constant in J/kg.K";
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta=u/d_crit "reduced density";
  Real tau=T_crit/T "inverse reduced temperature";

  SpecificEnthalpy h_of_u;

algorithm
  h_of_u := T*R*(tau*(EoS.f_it(
                           delta=delta, tau=tau) + EoS.f_rt(
                                                        delta=delta, tau=tau)) + (1+delta*EoS.f_rd(
                                                                                               delta=delta, tau=tau)));
  // return the RESidual
  y := h - h_of_u;
end setState_ThX_RES;
