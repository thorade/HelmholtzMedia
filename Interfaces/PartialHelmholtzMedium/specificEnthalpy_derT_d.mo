within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d "returns enthalpy derivative (dh/dT)@d=const"
  input ThermodynamicState state;
  output DerEnthalpyByTemperature dhTd;

import HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS.*;

protected
  HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    f := setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    dhTd := EoS.dhTd(f);
  elseif (state.phase == 2) then
    sat:=setSat_T(T=state.T);
    dhTd := 5000000000000;
  end if;
end specificEnthalpy_derT_d;
