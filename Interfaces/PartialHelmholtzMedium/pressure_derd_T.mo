within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function pressure_derd_T "returns pressure derivative (dp/dd)@T=const"
  input ThermodynamicState state;
  output DerPressureByDensity dpdT;

protected
  EoS.HelmholtzDerivs f;

algorithm
  if (state.phase == 1) then
    f:=EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    dpdT := EoS.dpdT(f);
  elseif (state.phase == 2) then
    dpdT := Modelica.Constants.small; // zero
  end if;
end pressure_derd_T;
