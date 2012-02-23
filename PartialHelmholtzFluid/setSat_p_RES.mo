within HelmholtzFluids.PartialHelmholtzFluid;
function setSat_p_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

  input AbsolutePressure p;

protected
  SaturationProperties sat;

algorithm
  sat := setSat_T(T=u);
  y := p - sat.psat;
end setSat_p_RES;
