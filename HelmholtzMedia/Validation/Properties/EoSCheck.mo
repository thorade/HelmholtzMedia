within HelmholtzMedia.Validation.Properties;
partial model EoSCheck "Single-phase and saturation properties against check values (e.g. computed with CoolProp)"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia
    constrainedby HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
  // single-phase states (T in K, d in kg/m3) and check values in SI units
  parameter Real T[:];
  parameter Real d[size(T, 1)];
  parameter Real p[size(T, 1)];
  parameter Real cv[size(T, 1)];
  parameter Real cp[size(T, 1)];
  parameter Real w[size(T, 1)];
  parameter Real h[size(T, 1)];
  parameter Real s[size(T, 1)];
  // saturation (T in K) and check values in SI units
  parameter Real Tsat[:];
  parameter Real psat[size(Tsat, 1)];
  parameter Real dl[size(Tsat, 1)];
  parameter Real dv[size(Tsat, 1)];
  parameter Real hl[size(Tsat, 1)];
  parameter Real hv[size(Tsat, 1)];
  parameter Real sv[size(Tsat, 1)];
  parameter Real tol=1e-5
    "relative tolerance (h and s relative to |value| + 1e5 J/kg and |value| + 1e2 J/(kg K))";
  Integer nFailed;

protected
  parameter String fileName = "eosCheck_" + Medium.mediumName + ".txt";
  Medium.ThermodynamicState state;
  Medium.SaturationProperties sat;

  function check
    input String name;
    input Real value;
    input Real ref;
    input Real tol "absolute tolerance";
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
    for i in 1:size(T, 1) loop
      state := Medium.setState_dTX(d=d[i], T=T[i], X=Medium.reference_X);
      Modelica.Utilities.Streams.print("T=" + String(T[i]) + " d=" + String(d[i]), fileName);
      nFailed := nFailed + check("p", Medium.pressure(state), p[i], tol*abs(p[i]), fileName);
      nFailed := nFailed + check("cv", Medium.specificHeatCapacityCv(state), cv[i], tol*abs(cv[i]), fileName);
      nFailed := nFailed + check("cp", Medium.specificHeatCapacityCp(state), cp[i], tol*abs(cp[i]), fileName);
      nFailed := nFailed + check("w", Medium.velocityOfSound(state), w[i], tol*abs(w[i]), fileName);
      nFailed := nFailed + check("h", Medium.specificEnthalpy(state), h[i], tol*(abs(h[i]) + 1e5), fileName);
      nFailed := nFailed + check("s", Medium.specificEntropy(state), s[i], tol*(abs(s[i]) + 1e2), fileName);
    end for;
    for i in 1:size(Tsat, 1) loop
      sat := Medium.setSat_T(T=Tsat[i]);
      Modelica.Utilities.Streams.print("Tsat=" + String(Tsat[i]), fileName);
      nFailed := nFailed + check("psat", sat.psat, psat[i], tol*abs(psat[i]), fileName);
      nFailed := nFailed + check("dl", sat.liq.d, dl[i], tol*abs(dl[i]), fileName);
      nFailed := nFailed + check("dv", sat.vap.d, dv[i], tol*abs(dv[i]), fileName);
      nFailed := nFailed + check("hl", sat.liq.h, hl[i], tol*(abs(hl[i]) + 1e5), fileName);
      nFailed := nFailed + check("hv", sat.vap.h, hv[i], tol*(abs(hv[i]) + 1e5), fileName);
      nFailed := nFailed + check("sv", sat.vap.s, sv[i], tol*(abs(sv[i]) + 1e2), fileName);
    end for;
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, getInstanceName() + ": " + String(nFailed) + " check values not met, see " + fileName);
  end when;

  annotation (Documentation(info="<html>
<p>
Base model: evaluates pressure, isochoric and isobaric heat capacity, speed of sound, enthalpy and entropy at single-phase states (T, d)
and the saturation properties at given temperatures, compares them with check values and writes the results to
<code>eosCheck_&lt;mediumName&gt;.txt</code>. Since the check values were computed with a different gas constant
(the one of the original equation of state) or are given with 7 significant digits, the tolerance is relative (default 1e-5).
</p>
</html>"));
end EoSCheck;
