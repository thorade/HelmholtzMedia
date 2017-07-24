within HelmholtzMedia.Examples.Validation;
model CriticalState_approach
  replaceable package Medium = HelmholtzFluids.Carbondioxide;

  Medium.ThermodynamicState state_sup;
  Medium.ThermodynamicState state_sub;

  Medium.EoS.HelmholtzDerivs f_sup = Medium.EoS.setHelmholtzDerivsThird(T=state_sup.T, d=state_sup.d, phase=state_sup.phase);
  Medium.EoS.HelmholtzDerivs f_sub = Medium.EoS.setHelmholtzDerivsThird(T=state_sub.T, d=state_sub.d, phase=state_sub.phase);

Modelica.Blocks.Sources.Ramp ramp(          duration=1,
    height=-1,
    offset=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));

protected
  constant Medium.Density dcrit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature Tcrit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.AbsolutePressure pcrit=Medium.fluidConstants[1].criticalPressure;
  constant Medium.SpecificEnthalpy HCRIT0=Medium.fluidConstants[1].HCRIT0;
  constant Medium.SpecificEntropy SCRIT0=Medium.fluidConstants[1].SCRIT0;

equation
  state_sup = Medium.setState_dT(d=dcrit+ramp.y, T=Tcrit, phase=1);
  state_sub = Medium.setState_dT(d=dcrit-ramp.y, T=Tcrit, phase=1);
end CriticalState_approach;
