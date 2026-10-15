within HelmholtzMedia.Validation.Properties;
model TransportPropylene "Propylene thermal conductivity against the check value of Assael et al. (2016), without critical enhancement"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Propylene;

  constant Medium.Temperature T_lambda = 350 "thermal conductivity check state";
  constant Medium.Density d_lambda = 385.0 "thermal conductivity check state";
  constant Real lambda_ref = 23.07 + 53.88
    "Assael et al., J. Phys. Chem. Ref. Data 45, 033104 (2016), section 3.2.4: dilute gas plus background contribution in mW/(m*K); the critical enhancement (4.52) is not included, since it needs the viscosity";
  constant Real tol_lambda = 2e-4 "relative tolerance (contributions given with 2 decimals)";
  Integer nFailed;

protected
  final constant String fileName = "transportPropylene.txt";
  Medium.ThermodynamicState state;
  Real lambda "thermal conductivity in mW/(m*K)";
  Real err;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    state := Medium.setState_dTX(d=d_lambda, T=T_lambda, X=Medium.reference_X);
    lambda := 1e3*Medium.thermalConductivity(state);
    err := lambda/lambda_ref - 1;
    if abs(err) > tol_lambda then
      nFailed := nFailed + 1;
    end if;
    Modelica.Utilities.Streams.print("T=" + String(T_lambda) + " d=" + String(d_lambda) + " lambda=" + String(lambda, significantDigits=10)
      + " ref=" + String(lambda_ref) + " rel_err=" + String(err, significantDigits=3), fileName);
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "TransportPropylene: " + String(nFailed) + " check value not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks the propylene thermal conductivity against the test point of Assael et al. (2016) at 350 K and 385 kg/m3.
Since HelmholtzMedia has no viscosity correlation for propylene, the critical enhancement is not computed, and the
check value is the sum of the dilute gas and background contributions (23.07 + 53.88 mW/(m K)).</p>
</html>"));
end TransportPropylene;
