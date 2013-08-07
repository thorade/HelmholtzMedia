within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function jouleThomsonCoefficient
  "returns Joule-Thomson-Coefficient (dT/dp)@h=const"
  input ThermodynamicState state;
  output DerTemperatureByPressure mu;

protected
  EoS.HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    mu := 1.0/(EoS.dpTd(f)-EoS.dpdT(f)*EoS.dhTd(f)/EoS.dhdT(f));
  elseif (state.phase == 2) then
    sat := setSat_T(T=state.T);
    mu := (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);
  end if;
end jouleThomsonCoefficient;
