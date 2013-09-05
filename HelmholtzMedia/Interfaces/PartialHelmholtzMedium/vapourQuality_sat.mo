within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function vapourQuality_sat "returns the vapour quality"

  input ThermodynamicState state;
  input SaturationProperties sat;
  output MassFraction x;

protected
  Temperature T_trip=fluidConstants[1].triplePointTemperature;
  Temperature T_crit=fluidConstants[1].criticalTemperature;

algorithm
  assert(state.T >= T_trip, "vapourQuality warning: Temperature is lower than triple-point temperature", level=AssertionLevel.warning);
  assert(state.T <= T_crit, "vapourQuality warning: Temperature is higher than critical temperature", level=AssertionLevel.warning);

  if state.d <= sat.vap.d then
    x := 1;
  elseif state.d >= sat.liq.d then
    x := 0;
  else
    x := (1.0/state.d - 1.0/sat.liq.d)/(1.0/sat.vap.d - 1.0/sat.liq.d);
  end if;

end vapourQuality_sat;
