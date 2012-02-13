within HelmholtzFluids.PartialHelmholtzFluid;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  constant Real[:,2] pressureSaturation "vapor pressure coefficients";
  constant Real[:,2] densityLiquid "saturated liquid density coefficients";
  constant Real[:,2] densityVapor "saturated vapor density coefficients";

end AncillaryCoefficients;
