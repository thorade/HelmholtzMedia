within HelmholtzMedia.Interfaces.Types;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  constant Temperature T_crit "critical temperature used by the ancillary equations (psat, dliq, dvap)";
  constant AbsolutePressure p_crit "critical pressure used by the ancillary equations (psat, dliq, dvap)";
  constant Density d_crit "critical density used by the ancillary equations (psat, dliq, dvap)";

  constant HelmholtzMedia.Interfaces.Choices.PressureSaturationModel pressureSaturationModel=HelmholtzMedia.Interfaces.Choices.PressureSaturationModel.PS5;
  constant Real[:,2] pressureSaturation = fill(0.0, 0, 2)
    "vapor pressure coefficients";

  constant HelmholtzMedia.Interfaces.Choices.DensityLiquidModel densityLiquidModel=HelmholtzMedia.Interfaces.Choices.DensityLiquidModel.DL1;
  constant Real[:,2] densityLiquid = fill(0.0, 0, 2)
    "saturated liquid density coefficients";

  constant HelmholtzMedia.Interfaces.Choices.DensityVaporModel densityVaporModel=HelmholtzMedia.Interfaces.Choices.DensityVaporModel.DV3;
  constant Real[:,2] densityVapor = fill(0.0, 0, 2)
    "saturated vapor density coefficients";

  constant HelmholtzMedia.Interfaces.Choices.PressureMeltingModel pressureMeltingModel=HelmholtzMedia.Interfaces.Choices.PressureMeltingModel.ML1;
  constant Temperature T_reducing=273.15;
  constant AbsolutePressure p_reducing=101325;
  constant Real[:,2] pressureMelting1 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting2 = fill(0.0, 0, 2)
    "melting pressure coefficients";
  constant Real[:,2] pressureMelting3 = fill(0.0, 0, 2)
    "melting pressure coefficients";
end AncillaryCoefficients;
