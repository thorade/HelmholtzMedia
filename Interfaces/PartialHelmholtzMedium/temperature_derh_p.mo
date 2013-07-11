within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_derh_p "returns temperature derivative (dT/dh)@p=const"
  extends Modelica.Icons.Function;
  input ThermodynamicState state "thermodynamic state record";
  output DerTemperatureByEnthalpy dThp "Temperature derivative w.r.t. enthalpy";

algorithm
  dThp := 1.0/specificHeatCapacityCp(state);
end temperature_derh_p;
