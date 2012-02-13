within HelmholtzFluids.PartialHelmholtzFluid;
function dewDensity_T_ANC
  "ancillary function: calculate density of saturated vapor for a given T"

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dvap;

protected
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real T_theta=1 - T/T_crit;
  Real tau=T_crit/T "inverse reduced temperature";

  Integer nDvap = size(ancillaryCoefficients.densityVapor,1);
  Real[nDvap] n = ancillaryCoefficients.densityVapor[:,1];
  Real[nDvap] theta = ancillaryCoefficients.densityVapor[:,2];

algorithm
    dvap := d_crit*exp(sum(n[i]*T_theta^theta[i] for i in 1:nDvap));  // DV3 Eric W. Lemmon
//  dvap := d_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:nDvap)); // DV6 Bücker & Wagner
end dewDensity_T_ANC;
