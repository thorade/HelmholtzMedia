within HelmholtzFluids.PartialHelmholtzFluid;
function ai_tau_tau "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_ideal_tau_tau
    "ideal part of dimensionless Helmholtz energy";

protected
  Integer nLog=size(helmholtzCoefficients.idealLog, 1);
  Integer nPower=size(helmholtzCoefficients.idealPower, 1);
  Integer nEinstein=size(helmholtzCoefficients.idealEinstein, 1);
//  Integer nCosh=size(helmholtzCoefficients.idealCosh, 1);
//  Integer nSinh=size(helmholtzCoefficients.idealSinh, 1);

  Real[nLog, 2] l=helmholtzCoefficients.idealLog;
  Real[nPower, 2] p=helmholtzCoefficients.idealPower;
  Real[nEinstein, 2] e=helmholtzCoefficients.idealEinstein;
//  Real[nCosh, 2] c=helmholtzCoefficients.idealCosh;
//  Real[nSinh, 2] s=helmholtzCoefficients.idealSinh;

algorithm
  alpha_ideal_tau_tau :=
      sum(l[i,1]*tau^(l[i,2]-2)*l[i,2]*(l[i,2]-1)/(tau^l[i,2]) - l[i,1]*tau^(l[i,2]-1)*l[i,2]^2/(tau^(l[i,2]+1)) for i in 1:nLog)
    + sum(p[i,1]*tau^(p[i,2]-2)*p[i,2]*(p[i,2]-1) for i in 1:nPower)
    - sum(e[i,1]*(-e[i,2])^2*exp(e[i,2]*tau)*(1 - exp(e[i,2]*tau))^(-2) for i in 1:nEinstein);
end ai_tau_tau;
