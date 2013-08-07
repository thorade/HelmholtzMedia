within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_pT_der "time derivative of specificEnthalpy_pT"
  extends Modelica.Icons.Function;
  input AbsolutePressure p;
  input Temperature T;
//input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  input ThermodynamicState state;
  input Real p_der "time derivative of pressure";
  input Real T_der "time derivative of temperature";
  output Real h_der "time derivative of specific Enthalpy";

algorithm
  h_der := p_der*isothermalThrottlingCoefficient(state=state)
         + T_der*specificHeatCapacityCp(state=state);

annotation (Inline=true);
end specificEnthalpy_pT_der;
