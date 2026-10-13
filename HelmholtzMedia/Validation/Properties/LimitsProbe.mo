within HelmholtzMedia.Validation.Properties;
model LimitsProbe "check that fluidLimits enclose the EoS states at the corners of the validity range"
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Butane;
  Medium.ThermodynamicState hot_dilute = Medium.setState_dTX(d=DMIN, T=TMAX, phase=1) "(TMAX, DMIN): highest entropy";
  Medium.ThermodynamicState cold_dense = Medium.setState_dTX(d=DMAX, T=TMIN, phase=1) "(TMIN, DMAX)";
  Medium.SaturationProperties sat "saturation at TMIN, solved in the algorithm";
  Medium.ThermodynamicState hot_pmax "(TMAX, PMAX), solved in the algorithm";
  Medium.ThermodynamicState cold_pmax "(TMIN, PMAX), solved in the algorithm";
  Medium.SaturationProperties sat_trip "saturation at the triple point temperature, solved in the algorithm";
protected
  constant Medium.Temperature TMIN = Medium.fluidLimits.TMIN;
  constant Medium.Temperature TMAX = Medium.fluidLimits.TMAX;
  constant Medium.Density DMIN = Medium.fluidLimits.DMIN;
  constant Medium.Density DMAX = Medium.fluidLimits.DMAX;
  constant Medium.AbsolutePressure PMAX = Medium.fluidLimits.PMAX;
  constant Medium.Temperature T_trip = Medium.fluidConstants[1].triplePointTemperature;
  constant Medium.AbsolutePressure p_trip = Medium.fluidConstants[1].triplePointPressure;
  final constant String fileName = "limits.txt";
algorithm
  when initial() then
    Modelica.Utilities.Files.remove(fileName);
    Modelica.Utilities.Streams.print("s_TMAX_DMIN=" + String(hot_dilute.s, significantDigits=6)
      + " h_TMAX_DMIN=" + String(hot_dilute.h, significantDigits=6)
      + " s_TMIN_DMAX=" + String(cold_dense.s, significantDigits=6)
      + " h_TMIN_DMAX=" + String(cold_dense.h, significantDigits=6)
      + " p_TMIN_DMAX=" + String(cold_dense.p, significantDigits=6), fileName);
    Modelica.Utilities.Streams.print("SMIN=" + String(Medium.fluidLimits.SMIN, significantDigits=6)
      + " SMAX=" + String(Medium.fluidLimits.SMAX, significantDigits=6)
      + " HMIN=" + String(Medium.fluidLimits.HMIN, significantDigits=6)
      + " HMAX=" + String(Medium.fluidLimits.HMAX, significantDigits=6)
      + " DMAX=" + String(DMAX, significantDigits=6)
      + " PMAX=" + String(PMAX, significantDigits=6), fileName);
    sat_trip := Medium.setSat_T(T=T_trip);
    Modelica.Utilities.Streams.print("p_trip_EoS=" + String(sat_trip.psat, significantDigits=6)
      + " p_trip=" + String(p_trip, significantDigits=6), fileName);
    // only checked if the triple point lies within the validity range of the EoS (T_trip >= TMIN),
    // otherwise the EoS is extrapolated there (e.g. Ethanol, Hexamethyldisiloxane)
    assert(T_trip < TMIN or abs(p_trip/sat_trip.psat - 1) < 0.05,
      "triplePointPressure differs by more than 5 % from the EoS, saturation pressure at T_trip is " + String(sat_trip.psat, significantDigits=6));
    sat := Medium.setSat_T(T=TMIN);
    Modelica.Utilities.Streams.print("h_liq_TMIN=" + String(sat.liq.h, significantDigits=6)
      + " s_liq_TMIN=" + String(sat.liq.s, significantDigits=6)
      + " d_liq_TMIN=" + String(sat.liq.d, significantDigits=6)
      + " d_vap_TMIN=" + String(sat.vap.d, significantDigits=6), fileName);
    // (T,p) corners last: setState_pT iterates, a failure here keeps the lines above
    hot_pmax := Medium.setState_pTX(p=PMAX, T=TMAX, phase=1);
    Modelica.Utilities.Streams.print("h_TMAX_PMAX=" + String(hot_pmax.h, significantDigits=6)
      + " s_TMAX_PMAX=" + String(hot_pmax.s, significantDigits=6)
      + " d_TMAX_PMAX=" + String(hot_pmax.d, significantDigits=6), fileName);
    cold_pmax := Medium.setState_pTX(p=PMAX, T=TMIN, phase=1);
    Modelica.Utilities.Streams.print("h_TMIN_PMAX=" + String(cold_pmax.h, significantDigits=6)
      + " s_TMIN_PMAX=" + String(cold_pmax.s, significantDigits=6)
      + " d_TMIN_PMAX=" + String(cold_pmax.d, significantDigits=6), fileName);

    assert(Medium.fluidLimits.SMAX >= hot_dilute.s,
      "SMAX is too low, should be at least s(TMAX,DMIN)=" + String(hot_dilute.s, significantDigits=6));
    assert(Medium.fluidLimits.SMIN <= min(sat.liq.s, cold_pmax.s),
      "SMIN is too high, should be at most " + String(min(sat.liq.s, cold_pmax.s), significantDigits=6));
    assert(Medium.fluidLimits.HMIN <= min(sat.liq.h, cold_pmax.h),
      "HMIN is too high, should be at most " + String(min(sat.liq.h, cold_pmax.h), significantDigits=6));
    assert(Medium.fluidLimits.HMAX >= max(hot_dilute.h, hot_pmax.h),
      "HMAX is too low, should be at least " + String(max(hot_dilute.h, hot_pmax.h), significantDigits=6));
    assert(DMAX >= max(sat.liq.d, cold_pmax.d),
      "DMAX is too low, should be at least " + String(max(sat.liq.d, cold_pmax.d), significantDigits=6));
  end when;

  annotation (Documentation(info="<html>
<p>Evaluates the EoS at the corners of the validity range given by <code>Medium.fluidLimits</code>
and checks that the limits SMIN, SMAX, HMIN, HMAX and DMAX enclose these states:</p>
<ul>
<li>SMAX &ge; s(TMAX, DMIN), the highest entropy (entropy increases with T and decreases with d)</li>
<li>SMIN &le; min(s of saturated liquid at TMIN, s(TMIN, PMAX))</li>
<li>HMIN &le; min(h of saturated liquid at TMIN, h(TMIN, PMAX))</li>
<li>HMAX &ge; max(h(TMAX, DMIN), h(TMAX, PMAX)), because above the inversion temperature h increases with p</li>
<li>DMAX &ge; max(d of saturated liquid at TMIN, d(TMIN, PMAX))</li>
<li>triplePointPressure agrees within 5 % with the saturation pressure of the EoS at triplePointTemperature (only if triplePointTemperature &ge; TMIN)</li>
</ul>
<p>Note that (TMIN, PMAX) lies in the solid region, so the EoS is extrapolated there and the resulting limits are conservative.
All values are also written to <code>limits.txt</code>, use <code>redeclare package Medium</code> to check other fluids.</p>
</html>"));
end LimitsProbe;
