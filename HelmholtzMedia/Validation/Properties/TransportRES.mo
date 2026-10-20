within HelmholtzMedia.Validation.Properties;
model TransportRES "Viscosity and thermal conductivity from residual entropy scaling against the sample values of Martinek et al. (2025) and Li et al. (2024)"
  extends Modelica.Icons.Example;
  package Propylene = HelmholtzMedia.HelmholtzFluids.Propylene;
  package MM = HelmholtzMedia.HelmholtzFluids.Hexamethyldisiloxane;

  constant Real tol_propylene = 1e-4
    "relative tolerance for propylene (the dilute gas viscosity of the publication comes from REFPROP, here from the RES polynomial)";
  constant Real tol_MM = 1e-4 "relative tolerance for the MM viscosity (dilute gas viscosity as for propylene)";
  constant Real tol_MM_lambda = 2e-3
    "relative tolerance for the MM thermal conductivity: Li et al. omit the critical enhancement for d > 2*d_crit, HelmholtzMedia does not (about 1.2e-3 here)";
  Integer nFailed;

protected
  final constant String fileName = "transportRES.txt";

  function check "compare, print and count"
    input String name;
    input Real T;
    input Real p;
    input Real value;
    input Real ref "published RES value";
    input Real experimental "experimental value";
    input Real tol;
    input String fileName;
    output Integer failed;
  protected
    Real err=value/ref - 1;
  algorithm
    failed := if abs(err) > tol then 1 else 0;
    Modelica.Utilities.Streams.print(name + " T=" + String(T) + " p=" + String(p) + " value=" + String(value, significantDigits=8)
      + " RES_ref=" + String(ref, significantDigits=8) + " rel_err=" + String(err, significantDigits=3)
      + " exp=" + String(experimental) + " dev_exp=" + String(value/experimental - 1, significantDigits=3)
      + (if failed > 0 then "  FAILED" else ""), fileName);
  end check;

algorithm
  when initial() then
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    // Martinek et al. (2025), supporting information, Samples_pure_fluids_output.txt: eta in micro Pa s
    nFailed := nFailed + check("Propylene eta", 88.7, 1.0,
      1e6*Propylene.dynamicViscosity(Propylene.setState_pTX(p=1.0, T=88.7, X=Propylene.reference_X)),
      13881.240, 14460.0, tol_propylene, fileName);
    nFailed := nFailed + check("MM eta", 293.15, 4304.0,
      1e6*MM.dynamicViscosity(MM.setState_pTX(p=4304.0, T=293.15, X=MM.reference_X)),
      522.092, 592.0, tol_MM, fileName);
    // Li et al. (2024), supporting information, Table_S5_SI.txt: lambda in W/(m K), including the critical enhancement
    nFailed := nFailed + check("MM lambda", 401.61, 4016e3,
      MM.thermalConductivity(MM.setState_pTX(p=4016e3, T=401.61, X=MM.reference_X)),
      0.08115, 0.08115, tol_MM_lambda, fileName);
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "TransportRES: " + String(nFailed) + " check values not met, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>
Compares the transport properties from residual entropy scaling (RES) with the sample values in the supporting information
of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity), which were computed with the equations of
state of REFPROP 10, which HelmholtzMedia also uses for both fluids. The experimental values given there are printed for information.
</p>
<ul>
<li>Propylene: tolerance 1e-4.</li>
<li>Hexamethyldisiloxane (MM): viscosity 1e-4; thermal conductivity 2e-3, because the code of Li et al. sets the critical
enhancement to zero for d &ge; 2&middot;d<sub>crit</sub> or T &gt; 1.4&middot;T<sub>crit</sub>, whereas HelmholtzMedia evaluates the
Olchowy-Sengers term without this cut-off (continuous, as in REFPROP). At the check state (654 kg/m3, d/d<sub>crit</sub> = 2.4)
the enhancement is 1.2e-4 W/(m K), 0.14&nbsp;% of the total.</li>
</ul>
</html>"));
end TransportRES;
