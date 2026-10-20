within HelmholtzMedia.Validation.Properties;
model TransportRESComparison
  "Compare viscosity and thermal conductivity from residual entropy scaling (RES) with the reference correlations of the same fluid on a (T, d) grid"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Propane
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium "fluid with its reference correlations";
  replaceable package MediumRES = HelmholtzMedia.HelmholtzFluids.Propane_RES
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium "same fluid with RES transport properties";

  parameter Integer nT_sub=8 "number of temperature intervals between 1.02*TMIN and 0.98*T_crit";
  parameter Integer nT_sup=8 "number of temperature intervals between 1.02*T_crit and TMAX";
  parameter Medium.AbsolutePressure p_max=min(70e6, Medium.fluidLimits.PMAX)
    "highest pressure of the compared states (the reference correlations are typically valid up to 50 to 100 MPa)";
  parameter Real d_red[:]={0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25}
    "densities relative to the critical density";
  parameter String fileName="transportRESComparison.txt";

  Integer n "number of compared states";
  Real AARD_eta "average absolute relative deviation of the RES viscosity";
  Real AARD_lambda "average absolute relative deviation of the RES thermal conductivity";
  Real max_eta "maximum absolute relative deviation of the RES viscosity";
  Real max_lambda "maximum absolute relative deviation of the RES thermal conductivity";

protected
  final constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  final constant Medium.Density d_crit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  final constant Medium.Temperature TMIN=Medium.fluidLimits.TMIN;
  final constant Medium.Temperature TMAX=Medium.fluidLimits.TMAX;
  final constant Medium.Density DMAX=Medium.fluidLimits.DMAX;
  Medium.Temperature T;
  Medium.Density d;
  Medium.SaturationProperties sat;
  Medium.ThermodynamicState state;
  MediumRES.ThermodynamicState stateRES;
  Real eta "viscosity in micro Pa s";
  Real eta_RES "viscosity in micro Pa s";
  Real lambda "thermal conductivity in mW/(m K)";
  Real lambda_RES "thermal conductivity in mW/(m K)";
  Real dev_eta;
  Real dev_lambda;
  Integer region "1 gas, 2 liquid, 3 supercritical, 4 near-critical";
  Boolean valid;

algorithm
  when initial() then
    Modelica.Utilities.Files.remove(fileName);
    Modelica.Utilities.Streams.print("T d p region eta eta_RES dev_eta lambda lambda_RES dev_lambda", fileName);
    n := 0;
    AARD_eta := 0;
    AARD_lambda := 0;
    max_eta := 0;
    max_lambda := 0;
    for i in 0:nT_sub + nT_sup + 1 loop
      T := if i <= nT_sub then 1.02*TMIN + (0.98*T_crit - 1.02*TMIN)*i/nT_sub
           else 1.02*T_crit + (TMAX - 1.02*T_crit)*(i - nT_sub - 1)/nT_sup;
      for j in 1:size(d_red, 1) loop
        d := d_red[j]*d_crit;
        valid := d < DMAX;
        region := 3;
        if valid and T < T_crit then
          sat := Medium.setSat_T(T=T);
          valid := d < sat.vap.d or d > sat.liq.d;
          region := if d < sat.vap.d then 1 else 2;
        end if;
        if valid then
          state := Medium.setState_dTX(d=d, T=T, X=Medium.reference_X, phase=1);
          valid := state.p > 0 and state.p <= p_max;
        end if;
        if valid then
          if abs(T/T_crit - 1) < 0.05 and abs(d/d_crit - 1) < 0.5 then
            region := 4;
          end if;
          stateRES := MediumRES.setState_dTX(d=d, T=T, X=MediumRES.reference_X, phase=1);
          eta := 1e6*Medium.dynamicViscosity(state);
          eta_RES := 1e6*MediumRES.dynamicViscosity(stateRES);
          lambda := 1e3*Medium.thermalConductivity(state);
          lambda_RES := 1e3*MediumRES.thermalConductivity(stateRES);
          dev_eta := eta_RES/eta - 1;
          dev_lambda := lambda_RES/lambda - 1;
          n := n + 1;
          AARD_eta := AARD_eta + abs(dev_eta);
          AARD_lambda := AARD_lambda + abs(dev_lambda);
          max_eta := max(max_eta, abs(dev_eta));
          max_lambda := max(max_lambda, abs(dev_lambda));
          Modelica.Utilities.Streams.print(String(T, significantDigits=6) + " " + String(d, significantDigits=6) + " "
            + String(state.p, significantDigits=6) + " " + String(region) + " "
            + String(eta, significantDigits=6) + " " + String(eta_RES, significantDigits=6) + " " + String(dev_eta, significantDigits=4) + " "
            + String(lambda, significantDigits=6) + " " + String(lambda_RES, significantDigits=6) + " " + String(dev_lambda, significantDigits=4), fileName);
        end if;
      end for;
    end for;
    AARD_eta := AARD_eta/max(n, 1);
    AARD_lambda := AARD_lambda/max(n, 1);
    Modelica.Utilities.Streams.print("n=" + String(n) + " AARD_eta=" + String(AARD_eta, significantDigits=3)
      + " max_eta=" + String(max_eta, significantDigits=3) + " AARD_lambda=" + String(AARD_lambda, significantDigits=3)
      + " max_lambda=" + String(max_lambda, significantDigits=3), fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>
Compares viscosity and thermal conductivity from residual entropy scaling (RES, package <code>MediumRES</code>, e.g.
<a href=\"modelica://HelmholtzMedia.HelmholtzFluids.Propane_RES\">Propane_RES</a>) with the reference correlations of the same
fluid (package <code>Medium</code>) on a grid of temperatures (1.02 TMIN to 0.98 T<sub>crit</sub> and 1.02 T<sub>crit</sub> to TMAX)
and densities (0.001 to 3.25 d<sub>crit</sub>), single-phase states with 0 &lt; p &le; p_max (default 70 MPa) only. Both packages use the same equation of state,
so the deviations are those of the transport models alone. The file <code>transportRESComparison.txt</code> lists every state
with its region (1 gas, 2 liquid, 3 supercritical, 4 near-critical: |T/T<sub>crit</sub> - 1| &lt; 0.05 and |d/d<sub>crit</sub> - 1| &lt; 0.5)
and the relative deviations; the last line gives the average and maximum absolute relative deviations.
There is no assert: this is a comparison, not a regression test. The script <code>compare_res_transport.py</code> runs it for
all fluids with a <code>_RES</code> variant and summarizes the deviations per region.
</p>
</html>"));
end TransportRESComparison;
