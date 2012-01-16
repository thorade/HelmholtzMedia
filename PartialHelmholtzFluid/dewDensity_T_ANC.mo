within HelmholtzFluids.PartialHelmholtzFluid;
function dewDensity_T_ANC
  "ancillary function: calculate density of saturated vapor for a given T"

  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.Density dvap;

protected
  Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real tau=T_crit/T "inverse reduced temperature";
  Real T_theta=1 - T/T_crit;
  Real[size(ancillaryCoefficients.n_dvap,1)] n=ancillaryCoefficients.n_dvap;
  Real[size(ancillaryCoefficients.theta_dvap,1)] theta=ancillaryCoefficients.theta_dvap;

algorithm
    dvap := d_crit*exp(tau*sum(n[i]*T_theta^theta[i] for i in 1:4));
end dewDensity_T_ANC;
