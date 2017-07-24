within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function temperature_pd_Waals
  input AbsolutePressure p;
  input Density d;
  output Temperature T;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R=fluidConstants[1].gasConstant/MM
    "specific gas constant";
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // van der Waals, as described by Span (2000)
  Real a = 27/64 * (R*R*T_crit*T_crit) /p_crit "correction for attraction";
  Real b = (R*T_crit)/(8*p_crit) "correction for volume";

algorithm
  // T := (p+a/v^2)*(v-b)/R;
  T := (p+a*d^2)*(1/d-b)/R;
  // Modelica.Utilities.Streams.print("van der Waals finished, T=" + String(T), "printlog.txt");

end temperature_pd_Waals;
