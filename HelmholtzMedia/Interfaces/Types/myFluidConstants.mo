within HelmholtzMedia.Interfaces.Types;
record myFluidConstants
  extends HelmholtzMedia.Interfaces.Types.TwoPhase.FluidConstants;
  Modelica.SIunits.MolarHeatCapacity gasConstant= Modelica.Constants.R;
end myFluidConstants;
