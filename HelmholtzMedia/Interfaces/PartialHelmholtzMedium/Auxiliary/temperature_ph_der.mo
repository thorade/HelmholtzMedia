within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function temperature_ph_der "time derivative of temperature_ph"

  input AbsolutePressure p;
  input SpecificEnthalpy h;
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real T_der "time derivative of temperature";

algorithm
  // in two-phase, T depends on p only: (dT/dh)@p = 0 (cp is infinite there)
  T_der := p_der*jouleThomsonCoefficient(state=state)
         + (if state.phase == 2 then 0 else h_der/specificHeatCapacityCp(state=state));

annotation (Inline=true);
end temperature_ph_der;
