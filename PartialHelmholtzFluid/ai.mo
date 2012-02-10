within HelmholtzFluids.PartialHelmholtzFluid;
function ai "ideal part of dimensionless Helmholtz energy"

  input Real delta;
  input Real tau;
  output Real alpha_ideal "ideal part of dimensionless Helmholtz energy";

protected
  Integer nIdeal=size(helmholtzCoefficients.ideal, 1);
  Real[nIdeal, 2] n=helmholtzCoefficients.ideal;

algorithm
  if delta>0 then
    alpha_ideal :=
      log(delta) + n[1,1]*log(tau) + n[2,1] + n[3,1]*tau
      + sum(n[i,1]*log(1 - exp(n[i,2]*tau)) for i in 4:nIdeal);
  else
    alpha_ideal := -Modelica.Constants.inf;
  end if;
end ai;
