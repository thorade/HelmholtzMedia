within HelmholtzFluids.PartialHelmholtzFluid;
record HelmholtzCoefficients
  "Coefficients for Helmholtz energy equations of state"

  //ideal gas part: substance specific coefficients
  constant Real[:,2] idealLog;
  constant Real[:,2] idealPower;
  constant Real[:,2] idealEinstein;
  // constant Real[:,2] idealCosh;
  // constant Real[:,2] idealSinh;

  //residual part: substance specific coefficients
  constant Real[:,4] residualPoly;
  constant Real[:,4] residualBwr;
  constant Real[:,12] residualGauss;

end HelmholtzCoefficients;
