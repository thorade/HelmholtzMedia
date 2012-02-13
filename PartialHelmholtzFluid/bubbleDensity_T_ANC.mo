within HelmholtzFluids.PartialHelmholtzFluid;
function bubbleDensity_T_ANC
  "ancillary function: calculate density of saturated liquid for a given T"

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dliq;

protected
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real T_theta=1 - T/T_crit;
  //

  Integer nDliq = size(ancillaryCoefficients.densityLiquid,1);
  Real[nDliq] n = ancillaryCoefficients.densityLiquid[:,1];
  Real[nDliq] theta = ancillaryCoefficients.densityLiquid[:,2];

algorithm
    dliq := d_crit*(1 + sum(n[i]*T_theta^theta[i] for i in 1:nDliq));
end bubbleDensity_T_ANC;
