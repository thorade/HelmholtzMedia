within HelmholtzMedia.Examples.Validation;
model CriticalState
  package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol;
  Medium.ThermodynamicState criticalState;

protected
  constant Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
  constant Medium.SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant Medium.Density d_crit=MM/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.SpecificEnthalpy HCRIT0=Medium.fluidConstants[1].HCRIT0;
  constant Medium.SpecificEntropy SCRIT0=Medium.fluidConstants[1].SCRIT0;
  constant Real eps=1e-6;

algorithm
  criticalState := Medium.setState_dTX(d=d_crit, T=T_crit, phase=1);
  assert(abs(criticalState.h-HCRIT0)<eps, "HCRIT0 is wrong, should be equal to criticalState.h=" + String(criticalState.h,significantDigits=15));
  assert(abs(criticalState.s-SCRIT0)<eps, "SCRIT0 is wrong, should be equal to criticalState.s=" + String(criticalState.s,significantDigits=15));

end CriticalState;
