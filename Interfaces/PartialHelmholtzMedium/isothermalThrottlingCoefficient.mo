within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient
  "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
//input HelmholtzDerivs is optional and will be used for single-phase only
  input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
output DerEnthalpyByPressure delta_T;

algorithm
  if (state.phase == 1) then
    // single phase definition as in Span(2000)
    delta_T := 1/state.d*(1-(1+f.delta*f.rd-f.delta*f.tau*f.rtd)/(1+2*f.delta*f.rd+f.delta^2*f.rdd));
  elseif (state.phase == 2) then
    delta_T := Modelica.Constants.inf; // divide by zero
  end if;
end isothermalThrottlingCoefficient;
