within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derd_T "returns pressure derivative (dp/dd)@T=const"
input ThermodynamicState state;
input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
output Types.DerPressureByDensity dpdT;

algorithm
  if (state.phase == 1) then
    dpdT := state.T*f.R*(1 + 2*f.delta*f.rd + f.delta^2*f.rdd);
  elseif (state.phase == 2) then
    dpdT := Modelica.Constants.small; // zero
  end if;
end pressure_derd_T;
