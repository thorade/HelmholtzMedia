within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function jouleThomsonCoefficient
  "returns Joule-Thomson-Coefficient (dT/dp)@h=const"
  input ThermodynamicState state;
  output DerTemperatureByPressure mu;

protected
  EoS.HelmholtzDerivs f;

algorithm
  if state.phase==1 then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    mu := 1.0/(EoS.dpTd(f)-EoS.dpdT(f)*EoS.dhTd(f)/EoS.dhdT(f));
  elseif state.phase==2 then
    mu := saturationTemperature_derp_sat(sat=setSat_T(T=state.T));
  end if;
end jouleThomsonCoefficient;
