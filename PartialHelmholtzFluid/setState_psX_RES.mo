within HelmholtzFluids.PartialHelmholtzFluid;
function setState_psX_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

  input AbsolutePressure p;
  input SpecificEntropy s;
  input FixedPhase phase=0;

algorithm
  // return the RESidual
  y := s - specificEntropy_pT(p=p,T=u,phase=phase);
end setState_psX_RES;
