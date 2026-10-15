within HelmholtzMedia.Validation.Properties;
model EoSAmmonia "Ammonia EoS against the check values of Gao et al. (2023), Table 7, and saturation values computed with CoolProp 8.0.0"
  extends EoSCheck(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia,
    T={200.0, 250.0, 320.0, 380.0, 405.0, 406.0, 275.0, 325.0, 375.0, 405.0, 520.0, 520.0, 520.0, 620.0, 620.0, 620.0},
    d=1e3*MM*{50.0, 42.0, 35.0, 26.0, 16.5, 13.696, 0.2, 0.95, 3.0, 10.0, 1.5, 13.696, 20.0, 5.5, 14.0, 45.0},
    p=1e6*{555.5631, 108.7263, 27.63275, 8.174128, 11.25287, 11.45246, 0.4270652, 2.093496, 6.044096, 11.2204, 5.96888, 35.83377, 56.3075, 24.38934, 58.52754, 1355.431},
    cv=1/MM*{56.15976, 50.25907, 46.95882, 47.82965, 77.13442, 103.0508, 32.60471, 40.6866, 48.9189, 77.15007, 36.52118, 43.66588, 43.81573, 41.29819, 43.50218, 50.78124},
    cp=1/MM*{74.22262, 72.04765, 78.92022, 125.1971, 4137.738, 52943.69, 45.39676, 65.25439, 111.3738, 2051.793, 49.38647, 98.36082, 87.32996, 61.81978, 75.46188, 64.36735},
    w={2940.432, 2049.168, 1370.446, 668.5456, 247.2696, 225.4458, 402.6253, 402.2104, 379.43, 286.2398, 538.8044, 555.4454, 741.1842, 588.0372, 692.5867, 3201.084},
    h=1e3/MM*{9.541211, 5.715934, 9.831787, 15.45317, 20.13541, 21.3041, 27.46383, 27.92681, 27.45404, 23.20689, 36.23166, 28.71997, 26.48937, 38.40927, 34.49909, 50.78325} .+ dh_IIR,
    s=1/MM*{-11.61817, 14.26145, 36.0605, 53.81385, 65.30611, 68.15535, 104.2956, 94.19375, 86.36105, 72.89756, 106.6103, 80.49762, 73.89537, 100.0682, 87.78122, 54.81426} .+ ds_IIR,
    Tsat={200, 240, 273.15, 300, 350, 400},
    psat={8609.776626, 102171.0302, 429247.7987, 1061121.502, 3865199.174, 10297199.4},
    dl={728.6693135, 681.4309219, 638.6382049, 600.1699235, 512.4162065, 344.0076036},
    dv={0.08867285906, 0.8969191288, 3.456010881, 8.24427316, 31.33643153, 130.8857835},
    hl={-126220.2624, 49471.08258, 200000, 326380.8388, 582657.1726, 942918.6633},
    hv={1351536.745, 1418647.623, 1461773.198, 1484432.155, 1478168.206, 1289109.072},
    sv={7003.620529, 6119.525895, 5619.341747, 5297.807952, 4769.385564, 3991.508427});

  // Table 7 is in molar units (density in mol/dm3, p in MPa, cv, cp, s in J/(mol K), h in kJ/mol)
  // with the reference state of Gao et al. (u = 0 and s = 0 for the saturated liquid at the triple point)
  constant Real MM = Medium.fluidConstants[1].molarMass;
  constant Real a1_Gao = -6.59406093943886 "a1 of Gao et al., Eq. (28)";
  constant Real a2_Gao = 5.60101151987913 "a2 of Gao et al., Eq. (28)";
  constant Real dh_IIR = Modelica.Constants.R/MM*Medium.helmholtzCoefficients.T_reducing*(Medium.helmholtzCoefficients.idealPower[2,1] - a2_Gao)
    "enthalpy offset IIR - Gao, in J/kg";
  constant Real ds_IIR = -Modelica.Constants.R/MM*(Medium.helmholtzCoefficients.idealPower[1,1] - a1_Gao)
    "entropy offset IIR - Gao, in J/(kg K)";
  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Single-phase states: the computer-verification values of Gao et al., J. Phys. Chem. Ref. Data 52, 013102 (2023), Table 7
(liquid, vapour and supercritical states, including states close to the critical point; the zero-density state is omitted),
converted to SI units and to the IIR reference state of this library. Gao et al. use R = 8.314462618 J/(mol K), the same as
Modelica.Constants.R, so the tolerance only reflects the 7 significant digits of the table.</p>
<p>Saturation: values computed with CoolProp 8.0.0 (same equation of state), IIR reference state.</p>
<p>Together with Validation.Derivatives.HelmholtzDerivsAssociating, this tests the associating terms of the residual Helmholtz energy.</p>
</html>"));
end EoSAmmonia;
