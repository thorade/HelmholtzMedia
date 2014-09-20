within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  constant PressureSaturationModel pressureSaturationModel=PressureSaturationModel.PS5;
  constant Real[:,2] pressureSaturation = fill(0.0, 0, 2)
    "vapor pressure coefficients";

  constant DensityLiquidModel densityLiquidModel=DensityLiquidModel.DL1;
  constant Real[:,2] densityLiquid = fill(0.0, 0, 2)
    "saturated liquid density coefficients";

  constant DensityVaporModel densityVaporModel=DensityVaporModel.DV3;
  constant Real[:,2] densityVapor = fill(0.0, 0, 2)
    "saturated vapor density coefficients";

  constant PressureMeltingModel pressureMeltingModel=PressureMeltingModel.ML1;
  constant Temperature T_reducing=273.15;
  constant AbsolutePressure p_reducing=101325;
  constant Real[:,2] pressureMelting1 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting2 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting3 = fill(0.0, 0, 2)
    "melting pressure coefficients";
end AncillaryCoefficients;
