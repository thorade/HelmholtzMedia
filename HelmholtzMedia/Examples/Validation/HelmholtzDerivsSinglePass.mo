within HelmholtzMedia.Examples.Validation;
model HelmholtzDerivsSinglePass "single-pass Helmholtz derivatives (f_i_derivs, f_r_derivs) against the separate functions"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Propane;
  parameter Integer nT=30;
  parameter Integer nd=30;
  parameter Real tol=1e-12 "relative tolerance for derivatives up to second order, scaled with the largest magnitude of each derivative on the grid";
  parameter Real tol3=2e-9
    "relative tolerance for the third derivatives: the expanded formulas in f_rttt, f_rddd lose up to about 1e-9 by cancellation";
  parameter String fileName="singlepass.txt";
  Integer nFailed;

protected
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Density d_crit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature T_min=max(Medium.fluidLimits.TMIN, 0.5*T_crit);
  constant Medium.Temperature T_max=min(Medium.fluidLimits.TMAX, 3*T_crit);
  constant Medium.Density d_min=1e-3*d_crit;
  constant Medium.Density d_max=min(Medium.fluidLimits.DMAX, 3*d_crit);
  constant String names[14]={"i","it","itt","ittt","r","rd","rdd","rt","rtt","rtd","rttt","rttd","rtdd","rddd"};
  Real delta;
  Real tau;
  Real new[14];
  Real old[14];
  Real maxAbs[14];
  Real maxDiff[14];

algorithm
  when initial() then
    maxAbs := zeros(14);
    maxDiff := zeros(14);
    for i in 1:nT loop
      tau := Medium.helmholtzCoefficients.T_reducing/(T_min*(T_max/T_min)^((i - 1)/(nT - 1)));
      for j in 1:nd loop
        delta := d_min*(d_max/d_min)^((j - 1)/(nd - 1))/Medium.helmholtzCoefficients.d_reducing;
        (new[1], new[2], new[3], new[4]) := Medium.EoS.f_i_derivs(delta=delta, tau=tau, order=3);
        (new[5], new[6], new[7], new[8], new[9], new[10], new[11], new[12], new[13], new[14]) :=
          Medium.EoS.f_r_derivs(delta=delta, tau=tau, order=3);
        old := {Medium.EoS.f_i(delta=delta, tau=tau), Medium.EoS.f_it(delta=delta, tau=tau),
                Medium.EoS.f_itt(delta=delta, tau=tau), Medium.EoS.f_ittt(delta=delta, tau=tau),
                Medium.EoS.f_r(delta=delta, tau=tau), Medium.EoS.f_rd(delta=delta, tau=tau),
                Medium.EoS.f_rdd(delta=delta, tau=tau), Medium.EoS.f_rt(delta=delta, tau=tau),
                Medium.EoS.f_rtt(delta=delta, tau=tau), Medium.EoS.f_rtd(delta=delta, tau=tau),
                Medium.EoS.f_rttt(delta=delta, tau=tau), Medium.EoS.f_rttd(delta=delta, tau=tau),
                Medium.EoS.f_rtdd(delta=delta, tau=tau), Medium.EoS.f_rddd(delta=delta, tau=tau)};
        for k in 1:14 loop
          maxAbs[k] := max(maxAbs[k], abs(old[k]));
          maxDiff[k] := max(maxDiff[k], abs(new[k] - old[k]));
        end for;
      end for;
    end for;
    nFailed := 0;
    Modelica.Utilities.Files.remove(fileName);
    for k in 1:14 loop
      if maxDiff[k] > (if k==4 or k>=11 then tol3 else tol)*max(maxAbs[k], 1) then
        nFailed := nFailed + 1;
      end if;
      Modelica.Utilities.Streams.print(names[k] + " maxAbs=" + String(maxAbs[k], significantDigits=6)
        + " maxDiff=" + String(maxDiff[k], significantDigits=3)
        + " rel=" + String(maxDiff[k]/max(maxAbs[k], 1), significantDigits=3), fileName);
    end for;
    Modelica.Utilities.Streams.print("nFailed=" + String(nFailed), fileName);
    assert(nFailed == 0, "HelmholtzDerivsSinglePass: " + String(nFailed) + " derivatives differ, see " + fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>Compares f_i_derivs and f_r_derivs (all derivatives in one pass) with the separate functions f_i, f_it, ..., f_rddd
on a grid of temperatures and densities. The difference of each derivative is scaled with its largest magnitude on the grid.</p>
<p>For the third derivatives the tolerance is larger: the expanded expressions in f_rttt and f_rddd suffer from cancellation
(e.g. Helium, Gaussian term with beta = -1357.66), which a high-precision evaluation (mpmath, 50 digits) showed to be the
less accurate of the two: absolute error 1.1e-6 of the separate function versus 9e-12 of f_r_derivs.</p>
</html>"));
end HelmholtzDerivsSinglePass;
