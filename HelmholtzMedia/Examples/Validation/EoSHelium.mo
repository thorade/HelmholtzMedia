within HelmholtzMedia.Examples.Validation;
model EoSHelium "Helium EoS against the check values of NIST IR 8474 (Tables 3 and 4); writes eosHelium.txt"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Helium;

  // Table 3: single-phase states, T in K, rho in mol/dm3, p in MPa, cv in J/(mol K), w in m/s
  constant Real T3[:] = {4, 4, 10, 10, 300, 300};
  constant Real rho3[:] = {40.0, 2.0, 50.0, 2.0, 25.0, 1.0};
  constant Real p3[:] = {1.593262, 0.0554523, 12.65519, 0.1588571, 85.769640, 2.524130};
  constant Real cv3[:] = {8.098737, 12.627957, 10.753076, 12.478387, 13.176055, 12.496256};
  constant Real w3[:] = {320.1490, 107.3812, 592.9440, 183.7793, 1349.3067, 1030.3609};
  // Table 4: saturation, T in K, psat in kPa, rho in mol/dm3, h in J/mol
  constant Real T4[:] = {2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.1, 5.15};
  constant Real psat4[:] = {5.3317, 8.3507, 12.376, 17.562, 24.061, 32.024, 41.599, 52.936, 66.186, 81.509, 99.076, 119.076, 141.732, 167.324, 196.235, 212.110, 220.461};
  constant Real rhol4[:] = {36.474, 36.321, 36.050, 35.697, 35.277, 34.796, 34.251, 33.638, 32.947, 32.164, 31.264, 30.209, 28.931, 27.296, 24.944, 23.109, 21.690};
  constant Real rhov4[:] = {0.30776, 0.45005, 0.62860, 0.84784, 1.1127, 1.4288, 1.8032, 2.2451, 2.7666, 3.3847, 4.1250, 5.0278, 6.1636, 7.6731, 9.9203, 11.705, 13.102};
  constant Real hl4[:] = {-27.1597, -24.9634, -22.9856, -20.9686, -18.8099, -16.4582, -13.8784, -11.0375, -7.8975, -4.4081, -0.4980, 3.9418, 9.0949, 15.3181, 23.5348, 29.4075, 33.6962};
  constant Real hv4[:] = {63.8799, 67.0948, 70.0799, 72.8144, 75.2766, 77.4416, 79.2798, 80.7528, 81.8099, 82.3803, 82.3597, 81.5857, 79.7820, 76.4226, 70.1954, 64.6840, 60.1771};
  // absolute tolerances: half a unit of the last tabulated digit plus 2e-6 relative,
  // which covers the different gas constant (8.314472 in the fit, Modelica.Constants.R here)
  constant Real tol_p3[:] = {3.69e-06, 1.61e-07, 3.03e-05, 3.68e-07, 0.000172, 5.55e-06};
  constant Real tol_cv3[:] = {1.67e-05, 2.58e-05, 2.2e-05, 2.55e-05, 2.69e-05, 2.55e-05};
  constant Real tol_w3[:] = {0.00069, 0.000265, 0.00124, 0.000418, 0.00275, 0.00211};
  constant Real tol_psat4[:] = {6.07e-05, 6.67e-05, 0.000525, 0.000535, 0.000548, 0.000564, 0.000583, 0.000606, 0.000632, 0.000663, 0.000698, 0.000738, 0.000783, 0.000835, 0.000892, 0.000924, 0.000941};
  constant Real tol_rhol4[:] = {0.000573, 0.000573, 0.000572, 0.000571, 0.000571, 0.00057, 0.000569, 0.000567, 0.000566, 0.000564, 0.000563, 0.00056, 0.000558, 0.000555, 0.00055, 0.000546, 0.000543};
  constant Real tol_rhov4[:] = {5.62e-06, 5.9e-06, 6.26e-06, 6.7e-06, 5.22e-05, 5.29e-05, 5.36e-05, 5.45e-05, 5.55e-05, 5.68e-05, 5.82e-05, 6.01e-05, 6.23e-05, 6.53e-05, 6.98e-05, 0.000523, 0.000526};
  constant Real tol_hl4[:] = {0.000104, 9.99e-05, 9.6e-05, 9.19e-05, 8.76e-05, 8.29e-05, 7.78e-05, 7.21e-05, 6.58e-05, 5.88e-05, 5.1e-05, 5.79e-05, 6.82e-05, 8.06e-05, 9.71e-05, 0.000109, 0.000117};
  constant Real tol_hv4[:] = {0.000178, 0.000184, 0.00019, 0.000196, 0.000201, 0.000205, 0.000209, 0.000212, 0.000214, 0.000215, 0.000215, 0.000213, 0.00021, 0.000203, 0.00019, 0.000179, 0.00017};
  Integer nFailed;

protected
  final constant String fileName = "eosHelium.txt";
  constant Real MM = Medium.fluidConstants[1].molarMass;
  Medium.ThermodynamicState state;
  Medium.SaturationProperties sat;
  Real x[5];

  function check
    input String name;
    input Real value;
    input Real ref;
    input Real tol;
    input String fileName;
    output Integer fail;
  algorithm
    fail := if abs(value - ref) > tol then 1 else 0;
    Modelica.Utilities.Streams.print("  " + name + "=" + String(value, significantDigits=10) + " ref=" + String(ref, significantDigits=10)
      + " dev=" + String(value - ref, significantDigits=3) + (if fail > 0 then "  FAILED" else ""), fileName);
  end check;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:size(T3, 1) loop
      state := Medium.setState_dTX(d=rho3[i]*1000*MM, T=T3[i], X=Medium.reference_X);
      x[1] := Medium.pressure(state)/1e6;
      x[2] := Medium.specificHeatCapacityCv(state)*MM;
      x[3] := Medium.velocityOfSound(state);
      Modelica.Utilities.Streams.print("T=" + String(T3[i]) + " rho=" + String(rho3[i]), fileName);
      nFailed := nFailed + check("p", x[1], p3[i], tol_p3[i], fileName);
      nFailed := nFailed + check("cv", x[2], cv3[i], tol_cv3[i], fileName);
      nFailed := nFailed + check("w", x[3], w3[i], tol_w3[i], fileName);
    end for;
    for i in 1:size(T4, 1) loop
      sat := Medium.setSat_T(T=T4[i]);
      x[1] := sat.psat/1e3;
      x[2] := sat.liq.d/(1000*MM);
      x[3] := sat.vap.d/(1000*MM);
      x[4] := sat.liq.h*MM;
      x[5] := sat.vap.h*MM;
      Modelica.Utilities.Streams.print("Tsat=" + String(T4[i]), fileName);
      nFailed := nFailed + check("psat", x[1], psat4[i], tol_psat4[i], fileName);
      nFailed := nFailed + check("rhol", x[2], rhol4[i], tol_rhol4[i], fileName);
      nFailed := nFailed + check("rhov", x[3], rhov4[i], tol_rhov4[i], fileName);
      nFailed := nFailed + check("hl", x[4], hl4[i], tol_hl4[i], fileName);
      nFailed := nFailed + check("hv", x[5], hv4[i], tol_hv4[i], fileName);
    end for;
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "EoSHelium: " + String(nFailed) + " check values not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks the Helium equation of state against the computer-verification values in Tables 3 and 4 of
Ortiz Vega et al., NIST IR 8474 (2023), DOI 10.6028/NIST.IR.8474.
The tolerance is half a unit of the last tabulated digit plus 2e-6 relative,
which allows for the gas constant (8.314472 J/(mol K) in the fit, Modelica.Constants.R in HelmholtzMedia).</p>
</html>"));
end EoSHelium;
