within HelmholtzMedia.Validation.Properties;
model TransportPropylene "Propylene thermal conductivity against the check value of Assael et al. (2016)"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Propylene;

  constant Medium.Temperature T_lambda = 350 "thermal conductivity check state";
  constant Medium.Density d_lambda = 385.0 "thermal conductivity check state";
  constant Real lambda_b_ref = 23.07 + 53.88
    "Assael et al., J. Phys. Chem. Ref. Data 45, 033104 (2016), section 3.2.4: dilute gas plus background contribution in mW/(m*K)";
  constant Real lambda_c_ref = 4.52 "critical enhancement in mW/(m*K), computed by Assael et al. with the REFPROP viscosity";
  constant Real tol_b = 2e-4 "relative tolerance of dilute gas plus background (contributions given with 2 decimals)";
  constant Real tol_total = 1e-2
    "relative tolerance of the total thermal conductivity (the critical enhancement uses the RES viscosity instead of the REFPROP viscosity)";
  Integer nFailed;

protected
  final constant String fileName = "transportPropylene.txt";
  Real lambda_b "dilute gas plus background in mW/(m*K)";
  Real lambda "total thermal conductivity in mW/(m*K)";
  Real err_b;
  Real err;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    lambda_b := 1e3*(Medium.Internal.Transport.thermalConductivity_dilute(d=d_lambda, T=T_lambda)
      + Medium.Internal.Transport.thermalConductivity_background(d=d_lambda, T=T_lambda));
    lambda := 1e3*Medium.thermalConductivity(Medium.setState_dTX(d=d_lambda, T=T_lambda, X=Medium.reference_X));
    err_b := lambda_b/lambda_b_ref - 1;
    err := lambda/(lambda_b_ref + lambda_c_ref) - 1;
    if abs(err_b) > tol_b then
      nFailed := nFailed + 1;
    end if;
    if abs(err) > tol_total then
      nFailed := nFailed + 1;
    end if;
    Modelica.Utilities.Streams.print("T=" + String(T_lambda) + " d=" + String(d_lambda)
      + " lambda_b=" + String(lambda_b, significantDigits=10) + " ref=" + String(lambda_b_ref) + " rel_err=" + String(err_b, significantDigits=3), fileName);
    Modelica.Utilities.Streams.print("T=" + String(T_lambda) + " d=" + String(d_lambda)
      + " lambda=" + String(lambda, significantDigits=10) + " ref=" + String(lambda_b_ref + lambda_c_ref) + " rel_err=" + String(err, significantDigits=3)
      + " lambda_c=" + String(lambda - lambda_b, significantDigits=4), fileName);
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "TransportPropylene: " + String(nFailed) + " check values not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks the propylene thermal conductivity against the test point of Assael et al. (2016) at 350 K and 385 kg/m3:
the dilute gas plus background contribution (23.07 + 53.88 mW/(m K)) to 2e-4, and the total including the critical
enhancement (4.52 mW/(m K)) to 1&nbsp;%, since the critical enhancement is computed with the viscosity from residual
entropy scaling instead of the REFPROP viscosity used by Assael et al.</p>
</html>"));
end TransportPropylene;
