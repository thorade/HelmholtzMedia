within HelmholtzMedia.Interfaces.Types;
record HelmholtzCoefficients
  "Coefficients for Helmholtz energy equations of state"

  //reducing parameters: tau=T_reducing/T, delta=d/d_reducing; usually, but not always, equal to the critical values
  constant Temperature T_reducing "reducing temperature of the EoS";
  constant Density d_reducing "reducing density of the EoS";

  //ideal gas part: substance specific coefficients
  constant Real[:,2] idealLog = fill(0.0, 0, 2);
  constant Real[:,2] idealPower = fill(0.0, 0, 2) "first two coeffs determine reference state";
  constant Real[:,2] idealEinstein = fill(0.0, 0, 2);
  constant Real[:,2] idealCosh = fill(0.0, 0, 2);
  constant Real[:,2] idealSinh = fill(0.0, 0, 2);

  //residual part: substance specific coefficients
  constant Real[:,3] residualPoly = fill(0.0, 0, 3);
  constant Real[:,4] residualBwr = fill(0.0, 0, 4);
  constant Real[:,9] residualGauss = fill(0.0, 0, 9);
  constant Real[:,12] residualNonAnalytical = fill(0.0, 0, 12);
  constant Real[:,8] residualAssociating = fill(0.0, 0, 8)
    "associating terms n*delta^d*tau^t*exp(eta*(delta-epsilon)^2 + 1/(beta*(tau-gamma)^2 + b)), columns <n, t, d, eta, beta, gamma, epsilon, b> (Gao et al. 2023, ammonia)";

  Boolean useLineSearch=false;
end HelmholtzCoefficients;
