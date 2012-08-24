within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_d_liq_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: d=u) and output y (Residual)

input Density d;

algorithm
  y := Ancillary.bubbleDensity_T(T=u) - d;
end saturationTemperature_d_liq_RES;
