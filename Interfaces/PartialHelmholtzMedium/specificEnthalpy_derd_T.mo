within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derd_T "returns enthalpy derivative (dh/dd)@T=const"
  input ThermodynamicState state;
//input HelmholtzDerivs is optional and will be used for single-phase only
  input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
//input sat is optional and will be used for two-phase only
  input SaturationProperties sat=setSat_T(T=state.T, phase=state.phase);
  output DerEnthalpyByDensity dhdT;

algorithm
  if (state.phase == 1) then
    dhdT := state.T*f.R/state.d*(0+f.tau*f.delta*f.rtd + f.delta*f.rd + f.delta^2*f.rdd);
  elseif (state.phase == 2) then
    // dhvT = (h"-h')/(v"-v')
    // dhdT = -1/d^2 * dhvT
    dhdT := -1/state.d^2*(sat.liq.h-sat.vap.h)/(1/sat.liq.d-1/sat.vap.d);
  end if;
end specificEnthalpy_derd_T;
