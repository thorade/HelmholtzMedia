within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_d_vap_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

input Density d;

algorithm
  y := Ancillary.dewDensity_T(T=u) - d;
end saturationTemperature_d_vap_RES;
