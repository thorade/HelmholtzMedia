within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
record HelmholtzCoefficients
  "Coefficients for Helmholtz energy equations of state"

  //ideal gas part: substance specific coefficients
  constant Real[:,2] idealLog = fill(0.0, 0, 2);
  constant Real[:,2] idealPower = fill(0.0, 0, 2);
  constant Real[:,2] idealEinstein = fill(0.0, 0, 2);
  constant Real[:,2] idealCosh = fill(0.0, 0, 2);
  constant Real[:,2] idealSinh = fill(0.0, 0, 2);

  //residual part: substance specific coefficients
  constant Real[:,4] residualPoly = fill(0.0, 0, 4);
  constant Real[:,4] residualBwr = fill(0.0, 0, 4);
  constant Real[:,9] residualGauss = fill(0.0, 0, 9);
  constant Real[:,12] residualNonAnalytical = fill(0.0, 0, 12);
end HelmholtzCoefficients;
