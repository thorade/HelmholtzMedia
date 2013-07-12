within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function temperature_ph_der "time derivative of temperature_ph"

  input AbsolutePressure p;
  input SpecificEnthalpy h;
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real T_der "time derivative of temperature";

algorithm
  T_der := p_der*jouleThomsonCoefficient(state=state)
         + h_der/specificHeatCapacityCp(state=state);

annotation (Inline=true);
end temperature_ph_der;
