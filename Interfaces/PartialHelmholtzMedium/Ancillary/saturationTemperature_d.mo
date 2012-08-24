within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_d
  "ancillary iterative function: calculate saturation temperature for a given density by iterating the ancillary function"

  input Modelica.SIunits.Density d;
  output Modelica.SIunits.Temperature T;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_trip=fluidConstants[1].triplePointTemperature;
  Temperature T_crit=fluidConstants[1].criticalTemperature;

  Real tolerance=1e-6 "tolerance for RES_d (in kg/m3)";

algorithm
  if (d-d_crit<tolerance) then
    // Modelica.Utilities.Streams.print("d<d_crit: vapour side", "printlog.txt");
    T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function saturationTemperature_d_vap_RES(d=d),
          u_min=T_trip,
          u_max=T_crit,
          tolerance=tolerance);
  elseif (d-d_crit>tolerance) then
    // Modelica.Utilities.Streams.print("d>d_crit: liquid side", "printlog.txt");
    T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function saturationTemperature_d_liq_RES(d=d),
          u_min=T_trip,
          u_max=T_crit,
          tolerance=tolerance);
  else
    // Modelica.Utilities.Streams.print("d=d_crit: return critical Temperature", "printlog.txt");
    T := T_crit;
  end if;

end saturationTemperature_d;
