within HelmholtzFluids.PartialHelmholtzFluid;
function dewDensity_T_ANC
  "ancillary function: calculate density of saturated vapor for a given T"
  import HelmholtzFluids.PartialHelmholtzFluid.Types.DensityVaporModel;

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dvap;

protected
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Real delta "reduced density";
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real T_theta;
  Real tau=T_crit/T "inverse reduced temperature";

  DensityVaporModel densityVaporModel=ancillaryCoefficients.densityVaporModel;
  Integer nDvap = size(ancillaryCoefficients.densityVapor,1);
  Real[nDvap] n = ancillaryCoefficients.densityVapor[:,1];
  Real[nDvap] theta = ancillaryCoefficients.densityVapor[:,2];

algorithm
    if (densityVaporModel == DensityVaporModel.DV1) then
      T_theta := (1 - T/T_crit); // odd
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := 1+delta; // DV1 or DV2
    elseif (densityVaporModel == DensityVaporModel.DV2) then
      T_theta := (1 - T/T_crit)^(1/3); // even
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := 1+delta; // DV1 or DV2

    elseif (densityVaporModel == DensityVaporModel.DV3) then
      T_theta := (1 - T/T_crit); // odd
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := exp(delta); // DV3 or DV4
    elseif (densityVaporModel == DensityVaporModel.DV4) then
      T_theta := (1 - T/T_crit)^(1/3); // even
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := exp(delta); // DV3 or DV4

    elseif (densityVaporModel == DensityVaporModel.DV5) then
      T_theta := (1 - T/T_crit); // odd
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := exp(tau*delta); // DV5 or DV6
    elseif (densityVaporModel == DensityVaporModel.DV6) then
      T_theta := (1 - T/T_crit)^(1/3); // even
      delta := sum(n[i]*T_theta^theta[i] for i in 1:nDvap);
      delta := exp(tau*delta); // DV5 or DV6
    end if;

    dvap := d_crit*delta;
end dewDensity_T_ANC;
