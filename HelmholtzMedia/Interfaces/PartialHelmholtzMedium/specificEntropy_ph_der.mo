within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEntropy_ph_der

  input AbsolutePressure p;
  input SpecificEnthalpy h;
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real h_der "time derivative of specific enthalpy";
  output Real s_der "time derivative of specific entropy";

algorithm
  s_der := p_der*(-1.0/(state.d*state.T))
         + h_der*(1.0/state.T);

annotation (
  Inline=true);
end specificEntropy_ph_der;
