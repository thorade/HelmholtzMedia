within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function pressure_dT_Waals
  input Density d;
  input Temperature T;
  output AbsolutePressure p;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=fluidConstants[1].gasConstant/MM
    "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // van der Waals, as described by Span (2000)
  Real a = 27/64 * (R_s*R_s*T_crit*T_crit) /p_crit "correction for attraction";
  Real b = (R_s*T_crit)/(8*p_crit) "correction for volume";

algorithm
  // p := R_s*T/(v-b) - a/v^2;
  p := R_s*T/(1/d-b) - a*d^2;
  // Modelica.Utilities.Streams.print("van der Waals finished, p=" + String(p), "printlog.txt");

end pressure_dT_Waals;
