within HelmholtzMedia.Examples.Validation;
model acentricFactor "validate acentric factor"
  package medium = HelmholtzFluids.Butane;

protected
  medium.Temperature T_crit = medium.fluidConstants[1].criticalTemperature;
  medium.AbsolutePressure p_crit = medium.fluidConstants[1].criticalPressure;
  Real omega1 = medium.fluidConstants[1].acentricFactor;
  Real omega2;

algorithm
  Modelica.Utilities.Streams.print("acentric Factor from parameter omega1=" + String(omega1), "printlog.txt");

  omega2 := -log10(HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary.saturationPressure_T(
                                             0.7*T_crit)/p_crit) - 1;
  Modelica.Utilities.Streams.print("acentric Factor from Wagner saturation pressure equation omega2=" + String(omega2), "printlog.txt");

annotation (experiment(NumberOfIntervals=1));
end acentricFactor;
