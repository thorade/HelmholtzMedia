within HelmholtzMedia.Examples.Validation;
model TransportR32 "R32 viscosity against the check values of Velliadou et al. (2022); writes transport properties to transportR32.txt"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.R32;

  constant Medium.Temperature T[:] = {300, 300, 300, 250, 400, 351.5, 360} "check states";
  constant Medium.Density d[:] = {1e-6, 10, 1100, 1150, 50, 424, 424} "check states";
  constant Real eta_ref[3] = {12.6170, 12.6333, 173.431}
    "Velliadou et al. (2022), Table 4: check values in muPa*s for the first three states";
  constant Real tol = 2e-5 "relative tolerance for the viscosity check values (6 significant digits)";
  Integer nFailed;

protected
  final constant String fileName = "transportR32.txt";
  Medium.ThermodynamicState state;
  Real eta "dynamic viscosity in muPa*s";
  Real lambda "thermal conductivity in mW/(m*K)";
  Real err;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:size(T, 1) loop
      state := Medium.setState_dTX(d=d[i], T=T[i], X=Medium.reference_X);
      eta := 1e6*Medium.dynamicViscosity(state);
      lambda := 1e3*Medium.thermalConductivity(state);
      if i <= size(eta_ref, 1) then
        err := eta/eta_ref[i] - 1;
        if abs(err) > tol then
          nFailed := nFailed + 1;
        end if;
      else
        err := 0;
      end if;
      Modelica.Utilities.Streams.print("T=" + String(T[i], significantDigits=15) + " d=" + String(d[i], significantDigits=15)
        + " eta=" + String(eta, significantDigits=12) + " lambda=" + String(lambda, significantDigits=12)
        + " eta_rel_err=" + String(err, significantDigits=6), fileName);
    end for;
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "TransportR32: " + String(nFailed) + " viscosity check values not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks the R32 viscosity against the computer-verification values of Velliadou et al. (2022)
at 300 K and 0, 10 and 1100 kg/m3 (the zero-density point is evaluated at 1e-6 kg/m3).
Viscosity and thermal conductivity at all states, including two near-critical ones,
are written to transportR32.txt for comparison with REFPROP.</p>
</html>"));
end TransportR32;
