within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d "returns enthalpy derivative (dh/dT)@d=const"
  input ThermodynamicState state;
  output DerEnthalpyByTemperature dhTd;

protected
  EoS.HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    f:=EoS.setHelmholtzDerivs(T=state.T, d=state.d, phase=state.phase);
    dhTd := f.R*(1 - f.tau^2*(f.itt+f.rtt) + f.delta*f.rd - f.tau*f.delta*f.rtd);
  elseif (state.phase == 2) then
    sat:=setSat_T(T=state.T);
    dhTd := 5000000000000;
  end if;
end specificEnthalpy_derT_d;
