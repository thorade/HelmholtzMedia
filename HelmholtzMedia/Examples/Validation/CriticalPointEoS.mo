within HelmholtzMedia.Examples.Validation;
model CriticalPointEoS "critical point of the EoS: (dp/dd)@T = (d2p/dd2)@T = 0"
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.Temperature T_c "critical temperature of the EoS";
  Medium.Density d_c "critical density of the EoS";
  Medium.AbsolutePressure p_c "critical pressure of the EoS";

protected
  constant Medium.Temperature T_crit = Medium.fluidConstants[1].criticalTemperature;
  constant Medium.AbsolutePressure p_crit = Medium.fluidConstants[1].criticalPressure;
  constant Medium.Density d_crit = Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  final constant String fileName = "criticalpoint.txt";

  function inflectionDensity "density where (d2p/dd2)@T = 0, by bisection"
    input Medium.Temperature T;
    output Medium.Density d;
  protected
    Medium.Density d_lo = 0.6*d_crit "gas-like side, d2p/dd2 < 0";
    Medium.Density d_hi = 1.5*d_crit "liquid-like side, d2p/dd2 > 0";
  algorithm
    for i in 1:100 loop
      d := (d_lo + d_hi)/2;
      if Medium.EoS.d2pd2T(Medium.EoS.setHelmholtzDerivsThird(d=d, T=T)) < 0 then
        d_lo := d;
      else
        d_hi := d;
      end if;
    end for;
  end inflectionDensity;

  function criticalPoint "temperature where the minimum of (dp/dd)@T, at the inflection density, is zero"
    output Medium.Temperature T;
    output Medium.Density d;
    output Medium.AbsolutePressure p;
  protected
    Medium.Temperature T_lo = 0.99*T_crit "dp/dd < 0 at the inflection density (van der Waals loop)";
    Medium.Temperature T_hi = 1.01*T_crit "dp/dd > 0 everywhere";
    Medium.EoS.HelmholtzDerivs f;
  algorithm
    assert(Medium.EoS.dpdT(Medium.EoS.setHelmholtzDerivsThird(d=inflectionDensity(T_lo), T=T_lo)) < 0
       and Medium.EoS.dpdT(Medium.EoS.setHelmholtzDerivsThird(d=inflectionDensity(T_hi), T=T_hi)) > 0,
       "CriticalPointEoS: 0.99*T_crit and 1.01*T_crit do not bracket the critical point of the EoS");
    for i in 1:100 loop
      T := (T_lo + T_hi)/2;
      d := inflectionDensity(T);
      f := Medium.EoS.setHelmholtzDerivsThird(d=d, T=T);
      if Medium.EoS.dpdT(f) < 0 then
        T_lo := T;
      else
        T_hi := T;
      end if;
    end for;
    p := Medium.EoS.p(f);
  end criticalPoint;

algorithm
  when initial() then
    (T_c, d_c, p_c) := criticalPoint();
    Modelica.Utilities.Files.remove(fileName);
    Modelica.Utilities.Streams.print("T_c=" + String(T_c, significantDigits=15) + " d_c=" + String(d_c, significantDigits=15)
      + " p_c=" + String(p_c, significantDigits=15), fileName);
    Modelica.Utilities.Streams.print("T_crit=" + String(T_crit, significantDigits=15) + " d_crit=" + String(d_crit, significantDigits=15)
      + " p_crit=" + String(p_crit, significantDigits=15), fileName);
    Modelica.Utilities.Streams.print("T_rel=" + String(T_c/T_crit-1) + " d_rel=" + String(d_c/d_crit-1) + " p_rel=" + String(p_c/p_crit-1), fileName);
    assert(abs(T_c/T_crit-1) < 1e-9 and abs(d_c/d_crit-1) < 1e-7 and abs(p_c/p_crit-1) < 1e-9,
      "critical constants differ from the critical point of the EoS, should be criticalTemperature=" + String(T_c, significantDigits=12)
      + ", criticalPressure=" + String(p_c, significantDigits=12) + ", critical density=" + String(d_c, significantDigits=12));
  end when;

  annotation (Documentation(info="<html>
<p>The critical constants in <code>fluidConstants</code> must be the critical point of the EoS itself,
where (dp/dd)@T = (d2p/dd2)@T = 0, because they are used for all region checks (two-phase or not).
They can differ from published critical values and from the reducing parameters of the EoS
(<code>helmholtzCoefficients.T_reducing</code>, <code>d_reducing</code>) and of the ancillary equations
(<code>ancillaryCoefficients.T_crit</code>, <code>p_crit</code>, <code>d_crit</code>).</p>
<p>This model determines the critical point of the EoS by two nested bisections: for a given temperature, the inflection
density with (d2p/dd2)@T = 0, and the temperature at which (dp/dd)@T at the inflection density is zero.
The results are written to <code>criticalpoint.txt</code>, and an assert checks that they agree with the critical constants
(relative tolerance 1e-9 for temperature and pressure, 1e-7 for density).</p>
</html>"));
end CriticalPointEoS;
