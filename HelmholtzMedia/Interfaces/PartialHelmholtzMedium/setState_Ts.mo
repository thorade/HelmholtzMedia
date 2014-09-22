within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_Ts "Return thermodynamic state as function of (T, s)"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input SpecificEntropy s "Entropy";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output ThermodynamicState state "thermodynamic state record";

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Temperature T_trip=fluidConstants[1].triplePointTemperature;

  EoS.HelmholtzDerivs f(T=T);
  SaturationProperties sat;
  MassFraction x "vapour quality";

  Density d_min=fluidLimits.DMIN;
  Density d_max=2*fluidLimits.DMAX; // extrapolation to higher densities should return reasonable values
  Density d_iter=d_crit;
  SpecificEntropy RES_s;
  SpecificEntropy RES_min;
  SpecificEntropy RES_max;
  DerEntropyByDensity dsdT "(ds/dd)@T=const";
  Real gamma(min=0.1,max=1) = 1 "convergence speed, default=1";
  constant Real tolerance=1e-9 "relative tolerance for RES_s (in J/kgK)";
  Integer iter = 0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print(" ", "printlog.txt");
  // Modelica.Utilities.Streams.print("setState_Ts: T=" + String(T) + " and s=" + String(s), "printlog.txt");
  state.phase := phase;

  if (state.phase == 2) then
    assert(T >= T_trip, "setState_Ts_error: pressure is lower than triple point pressure");
    assert(T <= T_crit, "setState_TsX_error: pressure is higher than critical pressure");
    sat := setSat_T(T=T);
    assert(s >= sat.liq.s, "setState_TsX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
    assert(s <= sat.vap.s, "setState_TsX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
  else
    if ((T < T_crit) and (T >= T_trip)) then
      // two-phase possible, check region
      if (T>0.95*T_crit) or (T<1.05*T_trip) then
        // Modelica.Utilities.Streams.print("close to critical or triple point, get saturation properties from EoS", "printlog.txt");
        sat := setSat_T(T=T);
      else
        // do a simple check first, quite often this is sufficient
        sat.liq.d := Ancillary.bubbleDensity_T(T=T);
        f := EoS.setHelmholtzDerivsFirst(T=T, d=sat.liq.d);
        sat.liq.s := EoS.s(f);

        sat.vap.d := Ancillary.dewDensity_T(T=T);
        f := EoS.setHelmholtzDerivsFirst(T=T, d=sat.vap.d);
        sat.vap.s := EoS.s(f);

        if ((s > sat.liq.s - abs(0.05*sat.liq.s)) and (s < sat.vap.s + abs(0.05*sat.vap.s))) then
          // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
          sat := setSat_T(T=T);
        end if;
      end if;

      if (s < sat.liq.s) then
        // Modelica.Utilities.Streams.print("single phase liquid: d is between dliq and rho_max", "printlog.txt");
        state.phase := 1;
        d_min := 0.98*sat.liq.d;
        d_iter := 1.02*sat.liq.d;
      elseif (s > sat.vap.s) then
        // Modelica.Utilities.Streams.print("single phase vapor: d is between 0 and sat.vap.d=" + String(sat.vap.d), "printlog.txt");
        state.phase := 1;
        d_max := sat.vap.d;
        d_iter := sat.vap.d/100;
      else
        // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
        state.phase := 2;
      end if;

    else
      // Modelica.Utilities.Streams.print("T>=T_crit or T<T_trip, only single phase possible, do not change dmin and dmax", "printlog.txt");
      state.phase := 1;
      d_iter := d_crit/100;
    end if;
  end if;

  // Modelica.Utilities.Streams.print("phase and region determination finished, d_min=" + String(d_min) + ", d_max=" + String(d_max) + " and d_iter=" + String(d_iter), "printlog.txt");

  if (state.phase == 2) then
    // force two-phase, SaturationProperties are already known
    state.T := T;
    state.s := s;
    state.p := sat.psat;
    x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
    state.d := 1.0/(1.0/sat.liq.d + x*(1.0/sat.vap.d - 1.0/sat.liq.d));
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
  else
    // force single-phase

    // Modelica.Utilities.Streams.print("initialize bisection", "printlog.txt");
    // min
    f := EoS.setHelmholtzDerivsFirst(T=T, d=d_min);
    RES_min := EoS.s(f) - s;
    // max
    f := EoS.setHelmholtzDerivsFirst(T=T, d=d_max);
    RES_max := EoS.s(f) - s;
    // iter
    f := EoS.setHelmholtzDerivsFirst(T=T, d=d_iter);
    RES_s := EoS.s(f) - s;

    //assert((RES_min*RES_max<0), "setState_Ts: d_min and d_max did not bracket the root, input was T=" + String(T) + " and s=" + String(s));
    // thighten the bounds
    // opposite sign brackets the root
    if (RES_s*RES_min<0) then
      d_max := d_iter;
      RES_max := RES_s;
    elseif (RES_s*RES_max<0) then
      d_min := d_iter;
      RES_min := RES_s;
    end if;

    while ((abs(RES_s/s) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;
      // gamma := iter/(iter+1);

      // calculate missing parts of f, then calculate gradient with respect to density
      f.rtd := EoS.f_rtd(delta=f.delta, tau=f.tau);
      dsdT := EoS.dsdT(f);

      // print for debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("d_iter=" + String(d_iter) + " and dsdT=" + String(dsdT), "printlog.txt");

      // calculate better d_iter
      d_iter := d_iter - gamma/dsdT*RES_s;

      // check bounds, if out of bounds use bisection
      if (d_iter<d_min) or (d_iter>d_max) then
        // Modelica.Utilities.Streams.print("d_iter out of bounds, fallback to bisection method, step=" + String(iter) + ", d_iter=" + String(d_iter), "printlog.txt");
        d_iter := (d_min+d_max)/2;
      end if;

      // calculate new RES_s
      f := EoS.setHelmholtzDerivsFirst(T=T, d=d_iter);
      RES_s := EoS.s(f) - s;

      // thighten the bounds
      // opposite sign brackets the root
      if (RES_s*RES_min<0) then
        d_max := d_iter;
        RES_max := RES_s;
      elseif (RES_s*RES_max<0) then
        d_min := d_iter;
        RES_min := RES_s;
      end if;

    end while;
    // Modelica.Utilities.Streams.print("setState_Ts total iteration steps " + String(iter), "printlog.txt");
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    assert(iter<iter_max, "setState_Ts did not converge, input was T=" + String(T) + " and s=" + String(s));

    state.T := T;
    state.s :=s;
    state.d := d_iter;
    state.p := EoS.p(f);
    state.h := EoS.h(f);
    state.u := EoS.u(f);
  end if;

end setState_Ts;
