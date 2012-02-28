within HelmholtzFluids.PartialHelmholtzFluid;
function jouleThomsonCoefficient "returns (dT/dp)@h=const"

input ThermodynamicState state;
output Types.JouleThomsonCoefficient mu;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Temperature T_trip=fluidConstants[1].triplePointTemperature;
  Real delta= state.d/d_crit "reduced density";
  Real tau= T_crit/state.T "inverse reduced temperature";

algorithm
  if (state.phase == 1) then
    // single phase definition as in Span(2000)
    mu := -(delta*ar_delta(delta=delta, tau=tau)
            + delta^2*ar_delta_delta(delta=delta, tau=tau)
            + delta*tau*ar_delta_tau(delta=delta, tau=tau))
            /
            ((1+delta*ar_delta(delta=delta, tau=tau))^2
            - tau^2*(ai_tau_tau(delta=delta, tau=tau)+ar_tau_tau(delta=delta, tau=tau))
            *(1+2*delta*ar_delta(delta=delta, tau=tau)+delta^2*ar_delta_delta(delta=delta, tau=tau)));
  elseif (state.phase == 2) then
    mu := Modelica.Constants.small; // zero? infinity?
  end if;
end jouleThomsonCoefficient;
