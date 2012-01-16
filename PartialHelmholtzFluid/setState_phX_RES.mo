within HelmholtzFluids.PartialHelmholtzFluid;
function setState_phX_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

  input AbsolutePressure p;
  input SpecificEnthalpy h;
  input FixedPhase phase=0;

algorithm
  // return the RESidual
  y := h - specificEnthalpy_pT(p=p,T=u,phase=phase);
end setState_phX_RES;
