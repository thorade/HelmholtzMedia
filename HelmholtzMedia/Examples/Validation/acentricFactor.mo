within HelmholtzMedia.Examples.Validation;
model acentricFactor "validate acentric factor"
  replaceable package Medium = HelmholtzFluids.Carbondioxide;

protected
  Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  Medium.AbsolutePressure p_crit=Medium.fluidConstants[1].criticalPressure;
  Real omega1=Medium.fluidConstants[1].acentricFactor;
  Real omega2;

algorithm
  omega2 := -log10(Medium.Ancillary.saturationPressure_T(0.7*T_crit)/p_crit) - 1;
  when terminal() then
    Modelica.Utilities.Streams.print("acentric Factor from parameter omega1=" + String(omega1), "printlog.txt");
    Modelica.Utilities.Streams.print("acentric Factor from Wagner saturation pressure equation omega2=" + String(omega2), "printlog.txt");
  end when;
end acentricFactor;
