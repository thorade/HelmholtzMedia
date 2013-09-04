within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function isentropicExpansionCoefficient "returns (dT/dp)@s=const"
  input ThermodynamicState state;
  output DerTemperatureByPressure delta_s;

protected
  EoS.HelmholtzDerivs f;

algorithm
  if (state.phase == 1) then
    f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
    delta_s := 1.0/(EoS.dpTd(f)-EoS.dpdT(f)*EoS.dsTd(f)/EoS.dsdT(f));
  elseif (state.phase == 2) then
    sat := setSat_T(T=state.T);
    delta_s := (1.0/sat.vap.d-1.0/sat.liq.d)/(sat.vap.s-sat.liq.s);
  end if;
end isentropicExpansionCoefficient;
