within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isothermalThrottlingCoefficient
  "returns isothermal throttling coefficient (dh/dp)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByPressure delta_T;

protected
  EoS.HelmholtzDerivs f;

algorithm
  if (state.phase == 1) then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
    delta_T := EoS.dhdT(f)/EoS.dpdT(f);
  elseif (state.phase == 2) then
    assert(false, "isothermalThrottlingCoefficient warning: in the two-phase region delta_T is infinite", level=AssertionLevel.warning);
    delta_T := Modelica.Constants.inf;
  end if;
end isothermalThrottlingCoefficient;
