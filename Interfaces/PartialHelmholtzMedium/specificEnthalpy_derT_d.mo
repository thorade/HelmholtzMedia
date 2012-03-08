within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d "returns enthalpy derivative (dh/dT)@d=const"
  input ThermodynamicState state;
//input HelmholtzDerivs is optional and will be used for single-phase only
  input HelmholtzDerivs f=setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
  output DerEnthalpyByTemperature dhTd;

algorithm
  if (state.phase == 1) then
    dhTd := f.R*(1 - f.tau^2*(f.itt+f.rtt) + f.delta*f.rd - f.tau*f.delta*f.rtd);
  elseif (state.phase == 2) then
    dhTd := 5000000000000;
  end if;
end specificEnthalpy_derT_d;
