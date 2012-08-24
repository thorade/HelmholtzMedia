within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function dewDensity_T
  "ancillary function: calculate density of saturated vapor for a given T"

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dvap;

protected
  DensityVaporModel densityVaporModel=ancillaryCoefficients.densityVaporModel;
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Real delta "reduced density";
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real T_theta;
  Real tau=T_crit/T "inverse reduced temperature";

  Integer nDvap=size(ancillaryCoefficients.densityVapor, 1);
  Real[nDvap] n=ancillaryCoefficients.densityVapor[:, 1];
  Real[nDvap] theta=ancillaryCoefficients.densityVapor[:, 2];

algorithm
  if (densityVaporModel == DensityVaporModel.DV1) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := 1 + delta; // DV1 or DV2
  elseif (densityVaporModel == DensityVaporModel.DV2) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := 1 + delta; // DV1 or DV2

  elseif (densityVaporModel == DensityVaporModel.DV3) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := exp(delta);   // DV3 or DV4
  elseif (densityVaporModel == DensityVaporModel.DV4) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := exp(delta);   // DV3 or DV4

  elseif (densityVaporModel == DensityVaporModel.DV5) then
    T_theta := max((1 - T/T_crit), Modelica.Constants.small); // odd
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := exp(tau*delta);   // DV5 or DV6
  elseif (densityVaporModel == DensityVaporModel.DV6) then
    T_theta := max((1 - T/T_crit)^(1/3), Modelica.Constants.small); // even
    delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
    delta := exp(tau*delta);   // DV5 or DV6
  end if;

  dvap := d_crit*delta;
end dewDensity_T;
