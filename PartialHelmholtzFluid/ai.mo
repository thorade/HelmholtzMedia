within HelmholtzFluids.PartialHelmholtzFluid;
function ai "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_ideal "ideal part of dimensionless Helmholtz energy";

protected
  Integer n_ideal=size(helmholtzCoefficients.n_ideal, 1);// not used yet, reorganize to Matrix structure
  Real[size(helmholtzCoefficients.n_ideal, 1)] n=helmholtzCoefficients.n_ideal;
  Real[size(helmholtzCoefficients.Theta, 1)] Theta=helmholtzCoefficients.Theta;

algorithm
  if delta>0 then
    alpha_ideal := log(delta) + n[1] + n[2]*tau + n[3]*log(tau)
      + sum(n[i]*log(1 - exp(-Theta[i]*tau)) for i in 4:7);
  else
    alpha_ideal := -Modelica.Constants.inf;
  end if;
end ai;
