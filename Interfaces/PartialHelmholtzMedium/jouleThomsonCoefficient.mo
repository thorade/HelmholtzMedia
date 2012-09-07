within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function jouleThomsonCoefficient
  "returns Joule-Thomson-Coefficient (dT/dp)@h=const"
  input ThermodynamicState state;
  output DerTemperatureByPressure mu;

protected
  EoS.HelmholtzDerivs f;
  SaturationProperties sat;

algorithm
  if (state.phase == 1) then
    f := EoS.setHelmholtzDerivsSecond(
                                T=state.T, d=state.d, phase=state.phase);
    // single phase definition as in Span(2000)
    mu := -(f.delta*f.rd + f.delta^2*f.rdd + f.delta*f.tau*f.rtd)/
           ((1+f.delta*f.rd)^2 - f.tau^2*(f.itt + f.rtt)*(1+2*f.delta*f.rd + f.delta^2*f.rdd));
  elseif (state.phase == 2) then
    sat := setSat_T(T=state.T);
    mu := saturationTemperature_derp(p=state.p, sat=sat);
  end if;
end jouleThomsonCoefficient;
