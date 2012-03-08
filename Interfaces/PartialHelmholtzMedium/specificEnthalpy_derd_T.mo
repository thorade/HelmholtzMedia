within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derd_T "returns enthalpy derivative (dh/dd)@T=const"
input ThermodynamicState state;
input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
output Types.DerEnthalpyByDensity dhdT;

protected
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    dhdT := state.T*f.R/state.d*(0+f.tau*f.delta*f.rtd + f.delta*f.rd + f.delta^2*f.rdd);
  elseif (state.phase == 2) then
    sat := setSat_T(T=state.T);
    dhdT := (sat.liq.h-sat.vap.h)/(sat.liq.d-sat.vap.d);
  end if;
end specificEnthalpy_derd_T;
