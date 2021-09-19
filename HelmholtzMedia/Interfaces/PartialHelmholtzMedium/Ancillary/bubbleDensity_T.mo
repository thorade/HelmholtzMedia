within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function bubbleDensity_T
  "ancillary function: calculate density of saturated liquid for a given T"

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dliq;

protected
  DensityLiquidModel densityLiquidModel=ancillaryCoefficients.densityLiquidModel;
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Real delta "reduced density";
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real T_theta;
  Real tau=T_crit/T "inverse reduced temperature";

  Integer nDliq=size(ancillaryCoefficients.densityLiquid, 1);
  Real[nDliq] n=ancillaryCoefficients.densityLiquid[:, 1];
  Real[nDliq] theta=ancillaryCoefficients.densityLiquid[:, 2];

algorithm
  if (densityLiquidModel == DensityLiquidModel.DL1) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := 1 + delta; // DL1 or DL2
  elseif (densityLiquidModel == DensityLiquidModel.DL2) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := 1 + delta; // DL1 or DL2

  elseif (densityLiquidModel == DensityLiquidModel.DL3) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := exp(delta);   // DL3 or DL4
  elseif (densityLiquidModel == DensityLiquidModel.DL4) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := exp(delta);   // DL3 or DL4

  elseif (densityLiquidModel == DensityLiquidModel.DL5) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := exp(tau*delta);   // DL5 or DL6
  elseif (densityLiquidModel == DensityLiquidModel.DL6) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta   := sum(n[i]*T_theta^theta[i] for i in 1:nDliq);
    delta   := exp(tau*delta);   // DL5 or DL6

  else
    assert(false, "bubbleDensity_d: this should not happen, check DensityLiquidModel");
  end if;

  dliq := d_crit*delta;

end bubbleDensity_T;
