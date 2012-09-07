within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derT_d "returns pressure derivative (dp/dT)@d=const"
  input ThermodynamicState state;
  output DerPressureByTemperature dpTd;

protected
  EoS.HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    dpTd := EoS.dpTd(f);
  elseif (state.phase == 2) then
    sat:=setSat_T(T=state.T);
    dpTd := saturationPressure_derT(T=state.T, sat=sat);
  end if;
end pressure_derT_d;
