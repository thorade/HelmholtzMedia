within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_pd "Return thermodynamic state as function of (p, d)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output ThermodynamicState state "thermodynamic state record";

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
    "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant Density d_reducing=helmholtzCoefficients.d_reducing "reducing density of the EoS";
  constant Temperature T_reducing=helmholtzCoefficients.T_reducing "reducing temperature of the EoS";
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
  constant Density dv_trip = Internal.Ancillary.dewDensity_T(T_trip);
  constant Density dl_trip = Internal.Ancillary.bubbleDensity_T(T_trip);

  Internal.EoS.HelmholtzDerivs f(d=d, delta = d/d_reducing);
  SaturationProperties sat;
  MassFraction x "vapour quality";

  Temperature T_min=0.98*fluidLimits.TMIN;
  Temperature T_max=2*max(fluidLimits.TMAX, p/(R_s*d));
  Temperature T_iter=T_crit;
  Temperature T_old "T_iter of the previous iteration step";
  AbsolutePressure RES_min;
  AbsolutePressure RES_max;
  AbsolutePressure RES_p;
  DerPressureByTemperature dpTd "(dp/dT)@d=const";
  Real gamma(min=0.1,max=1) = 1 "convergence speed, default=1";
  constant Real tolerance=1e-9 "relative tolerance for RES_p";
  Integer iter = 0;
  constant Integer iter_max = 200;
  Boolean converged = false "true once residual or Newton step is within tolerance";

algorithm
  // Modelica.Utilities.Streams.print(" ", "printlog.txt");
  // Modelica.Utilities.Streams.print("setState_pdX: p=" + String(p) + " and d=" + String(d), "printlog.txt");
  state.phase := phase;

  if state.phase==2 then
    assert(p <= p_crit, "setState_pdX_error: pressure is higher than critical pressure");
    sat := setSat_p(p=p);
    assert(d <= sat.liq.d, "setState_pdX_error: density is higher than saturated liquid density: this is single phase liquid");
    assert(d >= sat.vap.d, "setState_pdX_error: density is lower than saturated vapor density: this is single phase vapor");
  else
    if (p < p_crit) and (p>=p_trip) then
      // two-phase possible, do a region check
      if (p>0.98*p_crit) or (p<300*p_trip) then
        // Modelica.Utilities.Streams.print("close to critical or triple point", "printlog.txt");
        // get saturation properties from EoS
        sat := setSat_p(p=p);
      else
        // do a simple check first, quite often this is sufficient
        sat.Tsat := Internal.Ancillary.saturationTemperature_p(p=p);
        sat.liq.d := Internal.Ancillary.bubbleDensity_T(T=sat.Tsat);
        sat.vap.d := Internal.Ancillary.dewDensity_T(T=sat.Tsat);
        // Modelica.Utilities.Streams.print("setState_pd: sat.Tsat=" + String(sat.Tsat) + " and sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d) + ", simple check only", "printlog.txt");
        if ((d < sat.liq.d + abs(0.05*sat.liq.d)) and (d > sat.vap.d - abs(0.05*sat.vap.d))) then
          // Modelica.Utilities.Streams.print("setState_pd: p = " + String(p) + "d = " + String(d) + ", two-phase state or close to it", "printlog.txt");
          // get saturation properties from EoS
          sat := setSat_p(p=p);
        end if;
      end if;

      // Modelica.Utilities.Streams.print("setState_pd: phase boundary from EoS: sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d), "printlog.txt");
      if (d > sat.liq.d) then
        // Modelica.Utilities.Streams.print("single-phase liquid region", "printlog.txt");
        state.phase := 1;
        T_min := 0.98*Internal.Ancillary.saturationTemperature_d(d=d); // look at isobars in T,d-Diagram !!
        T_max := 1.02*sat.Tsat;
        T_iter := 1.1*T_min;
      elseif (d < sat.vap.d) then
        // Modelica.Utilities.Streams.print("single-phase vapour region", "printlog.txt");
        state.phase := 1;
        T_min := 0.99*sat.Tsat;
        T_iter:= Internal.Ancillary.temperature_pd_Waals(p=p, d=d);
      elseif (sat.liq.d > sat.vap.d) then
        // Modelica.Utilities.Streams.print("two-phase region, all properties can be calculated from sat record", "printlog.txt");
        state.phase := 2;
      else
        // liquid and vapour density coincide (at the critical point within roundoff): single phase, no vapour quality defined
        state.phase := 1;
        T_iter := sat.Tsat;
      end if;

    elseif (p < p_trip) then
      state.phase := 1;
      // very low pressure, behaves like an ideal gas
      T_iter := p/(d*R_s);

    elseif (p >= p_crit) then
      state.phase := 1;
      if (d>d_crit) then
        // Modelica.Utilities.Streams.print("p>p_crit and d>d_crit, single-phase super-critical liquid-like region", "printlog.txt");
        T_min := 0.98*Internal.Ancillary.saturationTemperature_d(d=d); // look at isobars in T,d-Diagram !!
        T_iter := 1.05*T_min;
      else
        // Modelica.Utilities.Streams.print("p>p_crit and d<=d_crit, single-phase super-critical vapour-like region", "printlog.txt");
        T_min := 0.98*T_crit;
        T_iter:= Internal.Ancillary.temperature_pd_Waals(p=p, d=d);
        T_max := 25*fluidLimits.TMAX;
      end if;
    else
      assert(false, "setState_pd: this should not happen, check p");
    end if;
  end if;
  // Modelica.Utilities.Streams.print("phase and region determination finshed, phase=" + String(state.phase) + ", T_min=" + String(T_min) + ", T_max=" + String(T_max) + ", T_iter=" + String(T_iter), "printlog.txt");

  // check bounds, the start values from van der Waals or the ancillary equations are not very accurate
  T_iter := max({T_iter, T_min, fluidLimits.TMIN});
  T_iter := min({T_iter, T_max, fluidLimits.TMAX});

  if state.phase==2 then
    // force two-phase, SaturationProperties are already known
    state.p := p;
    state.d := d;
    x := (1.0/d - 1.0/sat.liq.d)/(1.0/sat.vap.d - 1.0/sat.liq.d);
    state.T := sat.Tsat;
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
  else
    // force single-phase

    // Modelica.Utilities.Streams.print("initialize bisection", "printlog.txt");
    // min
    f := Internal.EoS.setHelmholtzDerivs(d=d, T=T_min, order=1, includeIdeal=false);
    RES_min := Internal.EoS.p(f) - p;
    // max
    f := Internal.EoS.setHelmholtzDerivs(d=d, T=T_max, order=1, includeIdeal=false);
    RES_max := Internal.EoS.p(f) - p;
    // iter
    f := Internal.EoS.setHelmholtzDerivs(d=d, T=T_iter, order=2, includeIdeal=false);
    RES_p := Internal.EoS.p(f) - p;

    assert((RES_min*RES_max<0), "setState_pd: T_min=" + String(T_min) + " and T_max=" + String(T_max) + " did not bracket the root, input was p=" + String(p) + " and d=" + String(d), level=AssertionLevel.warning);
    // thighten the bounds
    // opposite sign brackets the root
    if (RES_p*RES_min<0) then
      T_max := T_iter;
      RES_max := RES_p;
    elseif (RES_p*RES_max<0) then
      T_min := T_iter;
      RES_min := RES_p;
    end if;

    while (iter<iter_max) and not converged loop
      iter := iter+1;
      T_old := T_iter;
      // gamma := iter/(iter+1);

      // calculate gradients with respect to temperature (f contains the second derivatives of the residual part)
      dpTd := Internal.EoS.dpTd(f);

      // print for debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("T_iter=" + String(T_iter) + " and dpdT=" + String(dpdT), "printlog.txt");

      // calculate better T_iter
      T_iter := T_iter - gamma/dpTd*RES_p;

      // check bounds, if out of bounds use bisection
      if (T_iter<0.95*T_min) or (T_iter>1.05*T_max) then
        // Modelica.Utilities.Streams.print("T_iter out of bounds, fallback to bisection method, step=" + String(iter) + ", T_iter=" + String(T_iter), "printlog.txt");
        T_iter := (T_min+T_max)/2;
      end if;

      // calculate new RES_p
      f := Internal.EoS.setHelmholtzDerivs(d=d, T=T_iter, order=2, includeIdeal=false);
      RES_p := Internal.EoS.p(f) - p;

      // thighten the bounds
      // opposite sign brackets the root
      if (RES_p*RES_min<0) then
        T_max := T_iter;
        RES_max := RES_p;
      elseif (RES_p*RES_max<0) then
        T_min := T_iter;
        RES_min := RES_p;
      end if;

      // convergence test: relative residual within tolerance, or temperature step at machine precision
      // (in compressed liquid at low p, abs(RES_p/p)<=tolerance is not reachable)
      converged := abs(RES_p/p)<=tolerance or abs(T_iter-T_old)<=1e-12*T_iter;
    end while;
    // Modelica.Utilities.Streams.print("setState_pd required " + String(iter) + " iterations to find T=" + String(T_iter,significantDigits=6) + " for input p=" + String(p) + " and d=" + String(d), "printlog.txt");
    assert(converged, "setState_pdX did not converge (RES_p=" + String(RES_p) + "), input was p=" + String(p) + " and d=" + String(d));

    state.p := p;
    state.d := d;
    state.T := T_iter;
    f := Internal.EoS.setHelmholtzDerivs(d=state.d, T=state.T, order=1);
    state.h := Internal.EoS.h(f);
    state.u := Internal.EoS.u(f);
    state.s := Internal.EoS.s(f);
  end if;

end setState_pd;
