within HelmholtzFluids.PartialHelmholtzFluid;
function pressure_derd_T "returns (dp/dd)@T=const"

input ThermodynamicState state;
output Types.DerPressureByDensity pddT;

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
    pddT := R*state.T*(1 + 2*delta*ar_delta(delta=delta, tau=tau) + delta^2*ar_delta_delta(delta=delta, tau=tau));
  elseif (state.phase == 2) then
    pddT := Modelica.Constants.small; // zero
  end if;
end pressure_derd_T;
