within HelmholtzFluids.PartialHelmholtzFluid;
function pressure_derT_d "returns (dp/dT)@d=const"

input ThermodynamicState state;
output Types.DerPressureByTemperature pdTd;

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
    pdTd := R*state.d*(1 + delta*ar_delta(delta=delta, tau=tau) - delta*tau*ar_delta_tau(delta=delta, tau=tau));
  elseif (state.phase == 2) then
    pdTd := Modelica.Constants.inf; // How to calculate that?
  end if;
end pressure_derT_d;
