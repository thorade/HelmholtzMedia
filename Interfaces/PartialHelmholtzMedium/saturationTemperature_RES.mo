within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationTemperature_RES "residual function"
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
  // inherits input u (here: T=u) and output y (Residual)

  input AbsolutePressure p;

algorithm
  y := p - saturationPressure(T=u);
end saturationTemperature_RES;
