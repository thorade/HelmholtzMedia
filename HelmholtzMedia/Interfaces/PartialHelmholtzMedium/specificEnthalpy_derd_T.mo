within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derd_T "returns enthalpy derivative (dh/dd)@T=const"
  input ThermodynamicState state;
  output DerEnthalpyByDensity dhdT;

protected
  EoS.HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if state.phase==1 then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    dhdT := EoS.dhdT(f);
  elseif state.phase==2 then
    // dhvT = (h"-h')/(v"-v')
    // dhdT = -1/d^2 * dhvT
    sat:=setSat_T(T=state.T);
    dhdT := -1/state.d^2*(sat.liq.h-sat.vap.h)/(1/sat.liq.d-1/sat.vap.d);
  end if;
end specificEnthalpy_derd_T;
