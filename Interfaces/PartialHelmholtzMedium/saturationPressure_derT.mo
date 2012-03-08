within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function saturationPressure_derT "returns (dp/dT)@sat"

input Temperature T;
output DerPressureByTemperature dpT;

input SaturationProperties sat=setSat_T(T=T) "optional input sat";
// speeds up computation, if sat state is already known

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Temperature T_trip=fluidConstants[1].triplePointTemperature;
  Real tau= T_crit/sat.Tsat "inverse reduced temperature";
  Real delta_liq = sat.liq.d/d_crit;
  Real delta_vap = sat.vap.d/d_crit;

algorithm
    // algorithm by Span(2000) eq. 3.78
    dpT := (sat.vap.d*sat.liq.d)/(sat.vap.d-sat.liq.d)*R*(
        log(sat.vap.d/sat.liq.d)
        +     (f_r( tau=tau, delta=delta_vap)-f_r( tau=tau,delta=delta_liq))
        - tau*(f_rt(tau=tau, delta=delta_vap)-f_rt(tau=tau,delta=delta_liq)));
end saturationPressure_derT;
