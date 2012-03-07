within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  constant Types.PressureSaturationModel pressureSaturationModel;
  constant Real[:,2] pressureSaturation = fill(0.0, 0, 2)
    "vapor pressure coefficients";

  constant Types.DensityLiquidModel densityLiquidModel;
  constant Real[:,2] densityLiquid = fill(0.0, 0, 2)
    "saturated liquid density coefficients";

  constant Types.DensityVaporModel densityVaporModel;
  constant Real[:,2] densityVapor = fill(0.0, 0, 2)
    "saturated vapor density coefficients";
end AncillaryCoefficients;
