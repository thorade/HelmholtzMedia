within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function f_itt "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real f_ideal_tau_tau "ideal part of dimensionless Helmholtz energy";

protected
  constant Integer nLog=size(helmholtzCoefficients.idealLog, 1);
  constant Integer nPower=size(helmholtzCoefficients.idealPower, 1);
  constant Integer nEinstein=size(helmholtzCoefficients.idealEinstein, 1);
  constant Integer nCosh=size(helmholtzCoefficients.idealCosh, 1);
  constant Integer nSinh=size(helmholtzCoefficients.idealSinh, 1);

  constant Real[nLog, 2] l=helmholtzCoefficients.idealLog;
  constant Real[nPower, 2] p=helmholtzCoefficients.idealPower;
  constant Real[nEinstein, 2] e=helmholtzCoefficients.idealEinstein;
  constant Real[nCosh, 2] c=helmholtzCoefficients.idealCosh;
  constant Real[nSinh, 2] s=helmholtzCoefficients.idealSinh;

algorithm
  f_ideal_tau_tau :=
      sum(-(l[i,1]*l[i,2])/tau^2 for i in 1:nLog)
    + sum(p[i,1]*tau^(p[i,2]-2)*p[i,2]*(p[i,2]-1) for i in 1:nPower)
    - sum(e[i,1]*(-e[i,2])^2*exp(e[i,2]*tau)*(1 - exp(e[i,2]*tau))^(-2) for i in 1:nEinstein)
    - sum(c[i,1]*c[i,2]^2/(cosh(c[i,2]*tau))^2 for i in 1:nCosh)
    - sum(s[i,1]*s[i,2]^2/(sinh(s[i,2]*tau))^2 for i in 1:nSinh);
end f_itt;
