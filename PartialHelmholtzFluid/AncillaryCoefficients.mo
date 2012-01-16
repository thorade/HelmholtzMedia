within HelmholtzFluids.PartialHelmholtzFluid;
record AncillaryCoefficients
  "Coefficients for anciallry equations (psat, dliq, dvap)"

  // vapor pressure coefficients
  constant Real[:] n_vapor;
  constant Real[:] theta_vapor;

  // saturated liquid density coefficients
  constant Real[:] n_dliq;
  constant Real[:] theta_dliq;

  // saturated vapor density coefficients
  constant Real[:] n_dvap;
  constant Real[:] theta_dvap;

end AncillaryCoefficients;
