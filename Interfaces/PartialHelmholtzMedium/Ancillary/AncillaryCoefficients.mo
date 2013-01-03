within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
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

  constant Types.PressureMeltingModel pressureMeltingModel;
  constant Temperature T_reducing=1;
  constant AbsolutePressure p_reducing=1;
  constant Real[:,2] pressureMelting1 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting2 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting3 = fill(0.0, 0, 2)
    "melting pressure coefficients";
end AncillaryCoefficients;
