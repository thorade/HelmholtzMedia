within HelmholtzFluids.PartialHelmholtzFluid;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  constant Types.PressureSaturationModel pressureSaturationModel;
  constant Real[:,2] pressureSaturation "vapor pressure coefficients";

  constant Types.DensityLiquidModel densityLiquidModel;
  constant Real[:,2] densityLiquid "saturated liquid density coefficients";

  constant Types.DensityVaporModel densityVaporModel;
  constant Real[:,2] densityVapor "saturated vapor density coefficients";

end AncillaryCoefficients;
