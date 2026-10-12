within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient
  "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByPressure delta_T;

protected
  EoS.HelmholtzDerivs f;

algorithm
  if state.phase==1 then
    f := EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=1, order=2);
    delta_T := EoS.dhdT(f)/EoS.dpdT(f);
  elseif state.phase==2 then
    assert(false, "isothermalThrottlingCoefficient warning: in the two-phase region delta_T is infinite", level=AssertionLevel.warning);
    delta_T := Modelica.Constants.inf;
  else
    assert(false, "isothermalThrottlingCoefficient error: phase is neither 1 nor 2");
  end if;
end isothermalThrottlingCoefficient;
