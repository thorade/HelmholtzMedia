within HelmholtzFluids.PartialHelmholtzFluid;
record HelmholtzCoefficients
  "Coefficients for Helmholtz energy equations of state"

  //ideal gas part: substance specific coefficients
  constant Real[:,2] ideal;

  //residual part: substance specific coefficients
  constant Real[:,4] residualPoly;
  constant Real[:,4] residualBwr;
  constant Real[:,12] residualGauss;

end HelmholtzCoefficients;
