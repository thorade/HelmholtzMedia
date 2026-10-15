within HelmholtzMedia.Validation.Properties;
model TransportAmmonia "Ammonia viscosity and thermal conductivity against the check values of Monogenidou et al. (2018)"
  extends Modelica.Icons.Example;
  package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia;

  constant Medium.Temperature T_eta[:] = {300, 300, 300} "viscosity check states";
  constant Medium.Density d_eta[:] = {1e-6, 8.0, 609.0} "viscosity check states (zero density evaluated at 1e-6 kg/m3)";
  constant Real eta_ref[:] = {10.1812, 9.9219, 133.3937}
    "Monogenidou et al., J. Phys. Chem. Ref. Data 47, 023102 (2018), section 3: viscosity in muPa*s";
  constant Real tol_eta = 2e-5 "relative tolerance for the viscosity (4 decimals)";
  constant Medium.Temperature T_lambda = 390 "thermal conductivity check state";
  constant Medium.Density d_lambda = 415.0 "thermal conductivity check state";
  constant Real lambda_ref = 264.129743
    "Monogenidou et al., J. Phys. Chem. Ref. Data 47, 043101 (2018): thermal conductivity in mW/(m*K), including the critical enhancement 9.409965";
  constant Real tol_lambda = 1e-4
    "relative tolerance for the thermal conductivity (the critical enhancement depends on the EoS derivatives and the viscosity)";
  Integer nFailed;

protected
  final constant String fileName = "transportAmmonia.txt";
  Medium.ThermodynamicState state;
  Real eta "dynamic viscosity in muPa*s";
  Real lambda "thermal conductivity in mW/(m*K)";
  Real err;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for i in 1:size(T_eta, 1) loop
      state := Medium.setState_dTX(d=d_eta[i], T=T_eta[i], X=Medium.reference_X);
      eta := 1e6*Medium.dynamicViscosity(state);
      err := eta/eta_ref[i] - 1;
      if abs(err) > tol_eta then
        nFailed := nFailed + 1;
      end if;
      Modelica.Utilities.Streams.print("T=" + String(T_eta[i]) + " d=" + String(d_eta[i]) + " eta=" + String(eta, significantDigits=10)
        + " ref=" + String(eta_ref[i]) + " rel_err=" + String(err, significantDigits=3), fileName);
    end for;
    state := Medium.setState_dTX(d=d_lambda, T=T_lambda, X=Medium.reference_X);
    lambda := 1e3*Medium.thermalConductivity(state);
    err := lambda/lambda_ref - 1;
    if abs(err) > tol_lambda then
      nFailed := nFailed + 1;
    end if;
    Modelica.Utilities.Streams.print("T=" + String(T_lambda) + " d=" + String(d_lambda) + " lambda=" + String(lambda, significantDigits=10)
      + " ref=" + String(lambda_ref) + " rel_err=" + String(err, significantDigits=3), fileName);
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "TransportAmmonia: " + String(nFailed) + " check values not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Checks the ammonia viscosity against the computer-verification values of Monogenidou et al. (2018) at 300 K
and 0, 8 and 609 kg/m3, and the thermal conductivity (including the critical enhancement) at 390 K and 415 kg/m3.
The results are written to transportAmmonia.txt.</p>
</html>"));
end TransportAmmonia;
