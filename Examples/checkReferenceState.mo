within HelmholtzMedia.Examples;
model checkReferenceState
  package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.SpecificEnthalpy h_ref;
  Medium.SpecificEntropy s_ref;

protected
  Medium.SaturationProperties sat;
  Medium.Temperature T_IIR = 273.15; // 0°C;
  Medium.Temperature T_ASHRAE = 233.15; // -40°C;
  Medium.AbsolutePressure p_NBP = 101325; // 1.01325 bar = 1 atm

algorithm
  sat := Medium.setSat_T(T=T_IIR);
  // sat := Medium.setSat_T(T=T_ASHRAE);
  // sat := Medium.setSat_p(p=p_NBP);

  h_ref := sat.liq.h;
  s_ref := sat.liq.s;

end checkReferenceState;
