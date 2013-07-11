within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_derp_h "returns temperature derivative (dT/dp)@h=const"
  extends Modelica.Icons.Function;
  input ThermodynamicState state "thermodynamic state record";
  output DerTemperatureByPressure dTph "Temperature derivative w.r.t. pressure";

algorithm
  dTph := jouleThomsonCoefficient(state);
end temperature_derp_h;
