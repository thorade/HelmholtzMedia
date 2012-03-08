within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derT_d "returns pressure derivative (dp/dT)@d=const"
input ThermodynamicState state;
input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
output Types.DerPressureByTemperature dpTd;

algorithm
  if (state.phase == 1) then
    dpTd := state.d*f.R*(1 + f.delta*f.rd - f.delta*f.tau*f.rtd);
  elseif (state.phase == 2) then
    dpTd := saturationPressure_derT(T=state.T);
  end if;
end pressure_derT_d;
