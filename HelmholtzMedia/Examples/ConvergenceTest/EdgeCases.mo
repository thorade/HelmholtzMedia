within HelmholtzMedia.Examples.ConvergenceTest;
model EdgeCases "check the iterative medium functions at the edges of the validity range"
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;

protected
  constant Medium.Temperature T_crit = Medium.fluidConstants[1].criticalTemperature;
  constant Medium.AbsolutePressure p_crit = Medium.fluidConstants[1].criticalPressure;
  constant Medium.Density d_crit = Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature TMIN = Medium.fluidLimits.TMIN;
  constant Medium.Temperature TMAX = Medium.fluidLimits.TMAX;
  constant Medium.AbsolutePressure PMAX = Medium.fluidLimits.PMAX;
  constant Medium.SpecificHeatCapacity R_s = Modelica.Constants.R/Medium.fluidConstants[1].molarMass;

  constant Real tol_p = 2e-9 "relative tolerance for pressures, the solvers use 1e-9";
  constant Real tol_g = 1e-8 "tolerance for the Gibbs energy difference, relative to R_s*T";
  constant Modelica.Units.SI.Pressure p_floor = 1e-5 "absolute floor for liquid pressures at low p (floating-point resolution of the liquid)";
  constant Modelica.Units.SI.TemperatureDifference tol_T = 1e-6 "tolerance for round-trip temperatures";
  constant Real tol_d = 1e-8 "relative tolerance for round-trip densities";
  constant Real tol_s = 1e-9 "tolerance for entropies, relative to max(abs(s), R_s)";
  final constant String fileName = "edgecases.txt";

  Medium.SaturationProperties sat_TMIN "saturation at TMIN";
  Medium.SaturationProperties sat_low "saturation at p_low";
  Medium.SaturationProperties sat_mod "saturation at p_mod";
  Medium.SaturationProperties sat_a "saturation at T_a (bisection)";
  Medium.SaturationProperties sat_b "saturation at T_b (bisection)";
  Medium.AbsolutePressure p_low "low pressure, slightly above the saturation pressure at TMIN";
  Medium.AbsolutePressure p_mod "moderate pressure";
  Medium.Temperature T_liq "subcooled liquid at p_mod";
  Medium.Temperature T_vap "superheated vapour at p_mod";
  Medium.ThermodynamicState st;
  Medium.Temperature T_a;
  Medium.Temperature T_b;
  Medium.Temperature T_s0 "temperature where the saturated liquid entropy is zero";
  Integer nFailed(start=0, fixed=true) "number of failed checks";

  impure function report "write a failed check to the simulation log and to the result file"
    input String msg;
  algorithm
    Modelica.Utilities.Streams.print("EdgeCases FAILED: " + msg);
    Modelica.Utilities.Streams.print("FAILED: " + msg, fileName);
  end report;

  function liquidOk "liquid pressure equals p within tol_p, or as density error within tol_p (as in the solvers), or within p_floor"
    input Medium.ThermodynamicState liq;
    input Medium.AbsolutePressure p;
    output Boolean ok;
  algorithm
    ok := abs(liq.p/p-1)<=tol_p or abs((liq.p-p)/Medium.pressure_derd_T(liq))<=tol_p*liq.d or abs(liq.p-p)<=p_floor;
  end liquidOk;

  impure function checkSat_p "setSat_p: liquid and vapour pressure equal p, equal Gibbs energy"
    input Medium.AbsolutePressure p;
    input String name;
    output Boolean ok = true;
  protected
    Medium.SaturationProperties sat = Medium.setSat_p(p=p);
    Real gl = sat.liq.h - sat.Tsat*sat.liq.s;
    Real gv = sat.vap.h - sat.Tsat*sat.vap.s;
  algorithm
    Modelica.Utilities.Streams.print(name + ": Tsat=" + String(sat.Tsat, significantDigits=12) + " pl/p-1=" + String(sat.liq.p/p-1)
      + " pv/p-1=" + String(sat.vap.p/p-1) + " gl-gv=" + String(gl-gv), fileName);
    if not liquidOk(sat.liq, p) or abs(sat.vap.p/p-1)>tol_p or abs(gl-gv)>tol_g*R_s*sat.Tsat then
      report(name + ": inconsistent saturation state, pl/p-1=" + String(sat.liq.p/p-1) + ", pv/p-1=" + String(sat.vap.p/p-1) + ", gl-gv=" + String(gl-gv));
      ok := false;
    end if;
  end checkSat_p;

  impure function checkSat_T "setSat_T: equal liquid and vapour pressure, equal Gibbs energy"
    input Medium.Temperature T;
    input String name;
    output Boolean ok = true;
  protected
    Medium.SaturationProperties sat = Medium.setSat_T(T=T);
    Real gl = sat.liq.h - T*sat.liq.s;
    Real gv = sat.vap.h - T*sat.vap.s;
  algorithm
    Modelica.Utilities.Streams.print(name + ": psat=" + String(sat.psat, significantDigits=12) + " pl/pv-1=" + String(sat.liq.p/sat.vap.p-1)
      + " gl-gv=" + String(gl-gv), fileName);
    if not liquidOk(sat.liq, sat.vap.p) or abs(gl-gv)>tol_g*R_s*T then
      report(name + ": inconsistent saturation state, pl/pv-1=" + String(sat.liq.p/sat.vap.p-1) + ", gl-gv=" + String(gl-gv));
      ok := false;
    end if;
  end checkSat_T;

  impure function checkRoundTrip "pT -> ph and pT -> ps must return the original T and d"
    input Medium.AbsolutePressure p;
    input Medium.Temperature T;
    input String name;
    output Boolean ok = true;
  protected
    Medium.ThermodynamicState st = Medium.setState_pT(p=p, T=T);
    Medium.ThermodynamicState st_h = Medium.setState_ph(p=p, h=st.h);
    Medium.ThermodynamicState st_s = Medium.setState_ps(p=p, s=st.s);
  algorithm
    Modelica.Utilities.Streams.print(name + ": ph T_err=" + String(st_h.T-T) + " d_err=" + String(st_h.d/st.d-1)
      + " ps T_err=" + String(st_s.T-T) + " d_err=" + String(st_s.d/st.d-1), fileName);
    if abs(st_h.T-T)>tol_T or abs(st_h.d/st.d-1)>tol_d then
      report(name + ": pT->ph round trip, T_err=" + String(st_h.T-T) + ", d_err=" + String(st_h.d/st.d-1));
      ok := false;
    end if;
    if abs(st_s.T-T)>tol_T or abs(st_s.d/st.d-1)>tol_d then
      report(name + ": pT->ps round trip, T_err=" + String(st_s.T-T) + ", d_err=" + String(st_s.d/st.d-1));
      ok := false;
    end if;
  end checkRoundTrip;

  impure function checkState_pd "setState_pd: the EoS at the returned T must reproduce p"
    input Medium.AbsolutePressure p;
    input Medium.Density d;
    input String name;
    output Boolean ok = true;
  protected
    Medium.ThermodynamicState st = Medium.setState_pd(p=p, d=d);
    Medium.ThermodynamicState st_check = Medium.setState_dT(d=d, T=st.T, phase=1);
    Medium.AbsolutePressure p_check = Medium.pressure(st_check);
  algorithm
    Modelica.Utilities.Streams.print(name + ": T=" + String(st.T, significantDigits=12) + " p/p_in-1=" + String(p_check/p-1), fileName);
    // for a compressed liquid at very low p, p is a tiny difference of large terms and only resolved to about p_floor;
    // as in liquidOk, the deviation is then also accepted when it corresponds to a density error within tol_p
    if not (abs(p_check/p-1)<=tol_p or abs((p_check-p)/Medium.pressure_derd_T(st_check))<=tol_p*d or abs(p_check-p)<=p_floor) then
      report(name + ": setState_pd, p/p_in-1=" + String(p_check/p-1));
      ok := false;
    end if;
  end checkState_pd;

  impure function checkState_Ts "setState_Ts: the EoS at the returned d must reproduce s"
    input Medium.Temperature T;
    input Medium.SpecificEntropy s;
    input String name;
    output Boolean ok = true;
  protected
    Medium.ThermodynamicState st = Medium.setState_Ts(T=T, s=s);
    Medium.SpecificEntropy s_check = Medium.specificEntropy(Medium.setState_dT(d=st.d, T=T, phase=1));
  algorithm
    Modelica.Utilities.Streams.print(name + ": d=" + String(st.d, significantDigits=12) + " s_err=" + String(s_check-s), fileName);
    if abs(s_check-s)>tol_s*max(abs(s), R_s) then
      report(name + ": setState_Ts, s_err=" + String(s_check-s));
      ok := false;
    end if;
  end checkState_Ts;

algorithm
  when initial() then
    Modelica.Utilities.Files.remove(fileName);

    // test points derived from the medium
    sat_TMIN := Medium.setSat_T(T=TMIN);
    p_low := 1.5*sat_TMIN.psat;
    sat_low := Medium.setSat_p(p=p_low);
    p_mod := max(p_crit/4, 10*p_low);
    sat_mod := Medium.setSat_p(p=p_mod);
    T_liq := max(TMIN, 0.98*sat_mod.Tsat);
    T_vap := min(TMAX, 1.02*sat_mod.Tsat);

    // saturation from p: near the triple point, moderate, and near the critical point
    nFailed := nFailed + (if checkSat_p(p_low, "setSat_p(1.5*psat(TMIN))") then 0 else 1);
    nFailed := nFailed + (if checkSat_p(p_mod, "setSat_p(p_mod)") then 0 else 1);
    nFailed := nFailed + (if checkSat_p(p_crit*(1-1e-3), "setSat_p(p_crit*(1-1e-3))") then 0 else 1);
    nFailed := nFailed + (if checkSat_p(p_crit*(1-1e-6), "setSat_p(p_crit*(1-1e-6))") then 0 else 1);
    nFailed := nFailed + (if checkSat_p(p_crit*(1-3e-7), "setSat_p(p_crit*(1-3e-7))") then 0 else 1);

    // saturation from T: at TMIN and approaching the critical point
    nFailed := nFailed + (if checkSat_T(TMIN, "setSat_T(TMIN)") then 0 else 1);
    nFailed := nFailed + (if checkSat_T(T_crit*(1-1e-3), "setSat_T(T_crit*(1-1e-3))") then 0 else 1);
    nFailed := nFailed + (if checkSat_T(T_crit*(1-1e-6), "setSat_T(T_crit*(1-1e-6))") then 0 else 1);
    nFailed := nFailed + (if checkSat_T(T_crit*(1-3e-7), "setSat_T(T_crit*(1-3e-7))") then 0 else 1);
    nFailed := nFailed + (if checkSat_T(T_crit*(1-1e-9), "setSat_T(T_crit*(1-1e-9))") then 0 else 1);

    // round trips pT -> ph and pT -> ps
    nFailed := nFailed + (if checkRoundTrip(p_low, TMIN, "liquid at p_low, TMIN") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(p_low, TMAX, "vapour at p_low, TMAX") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(1e-3*p_low, TMAX, "vapour at 1e-3*p_low, TMAX") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(p_mod, T_liq, "liquid at p_mod") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(p_mod, T_vap, "vapour at p_mod") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(1.02*p_crit, 1.002*T_crit, "near-critical supercritical") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(PMAX, T_crit, "PMAX, T_crit") then 0 else 1);
    nFailed := nFailed + (if checkRoundTrip(PMAX, TMAX, "PMAX, TMAX") then 0 else 1);

    // setState_pd: compressed liquid at low pressure, vapour, supercritical
    nFailed := nFailed + (if checkState_pd(p_low, 1.0001*sat_low.liq.d, "liquid at p_low") then 0 else 1);
    nFailed := nFailed + (if checkState_pd(p_low, 0.5*sat_low.vap.d, "vapour at p_low") then 0 else 1);
    nFailed := nFailed + (if checkState_pd(1.02*p_crit, d_crit, "supercritical at d_crit") then 0 else 1);

    // setState_Ts: subcooled liquid and superheated vapour at p_mod
    st := Medium.setState_pT(p=p_mod, T=T_liq);
    nFailed := nFailed + (if checkState_Ts(T_liq, st.s, "liquid at p_mod") then 0 else 1);
    st := Medium.setState_pT(p=p_mod, T=T_vap);
    nFailed := nFailed + (if checkState_Ts(T_vap, st.s, "vapour at p_mod") then 0 else 1);

    // setState_Ts with s=0: depending on the reference state, the saturated liquid entropy may cross zero;
    // find that temperature by bisection and test s=0 in the compressed liquid slightly above it
    T_a := TMIN;
    T_b := T_crit*(1-1e-3);
    sat_a := Medium.setSat_T(T=T_a);
    sat_b := Medium.setSat_T(T=T_b);
    if (sat_a.liq.s < 0) and (sat_b.liq.s > 0) then
      for i in 1:60 loop
        T_s0 := (T_a+T_b)/2;
        sat_a := Medium.setSat_T(T=T_s0);
        if (sat_a.liq.s < 0) then
          T_a := T_s0;
        else
          T_b := T_s0;
        end if;
      end for;
      nFailed := nFailed + (if checkState_Ts(T_s0 + 0.02*(T_crit-T_s0), 0, "liquid with s=0") then 0 else 1);
    else
      Modelica.Utilities.Streams.print("liquid with s=0: not tested, saturated liquid entropy does not cross zero", fileName);
    end if;

    assert(nFailed == 0, String(nFailed) + " edge case check(s) failed for " + Medium.fluidConstants[1].chemicalFormula
      + ", see the lines starting with EdgeCases FAILED in the simulation log, or " + fileName);
  end when;

  annotation (Documentation(info="<html>
<p>Calls the iterative medium functions at the edges of the validity range, where convergence problems
have occurred in the past: near the triple point, near the critical point, in the compressed liquid at very low
pressure, at PMAX, and with s=0 (depending on the reference state). All test points are derived from
<code>Medium.fluidConstants</code> and <code>Medium.fluidLimits</code>, so the model can be used for every
fluid via <code>redeclare package Medium</code>.</p>
<p>Checked are the consistency of saturation states (equal pressure and Gibbs energy of liquid and vapour),
pT&rarr;ph and pT&rarr;ps round trips, and that setState_pd and setState_Ts reproduce their inputs.
Each failed check is reported in the simulation log, the number of failures is reported in one final assert; all results are written to <code>edgecases.txt</code>.
A function that does not converge at all stops the simulation with its own error message.</p>
<p>For the almost incompressible liquid, a pressure deviation is also accepted if it corresponds to a density deviation
within the tolerance, as in the solvers; at very low pressure, pressures are only resolvable to about 1e-7 Pa,
therefore an absolute floor <code>p_floor</code> is accepted there.</p>
</html>"));
end EdgeCases;
