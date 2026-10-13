within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificEnthalpy_ps_der "time derivative of specificEnthalpy_ps"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEntropy s "Specific entropy";
  input ThermodynamicState state "state at (p, s)";
  input Real p_der "time derivative of p";
  input Real s_der "time derivative of s";
  output Real h_der "time derivative of h";
algorithm
  h_der := p_der/state.d + s_der*state.T;
annotation (Inline=true);
end specificEnthalpy_ps_der;
