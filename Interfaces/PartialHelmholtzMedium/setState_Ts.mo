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
  Real delta "reduced density";
  Real tau(unit="1")=T_crit/T "inverse reduced temperature";
  EoS.HelmholtzDerivs f;

  SaturationProperties sat;
  MassFraction x "vapour quality";

  Density d_min=fluidLimits.DMIN;
  Density d_max=fluidLimits.DMAX;
  Density d_med;
  Density d_iter=d_crit;
  Real RES_s;
  Real RES_min;
  Real RES_max;
  Real RES_med;
  Real dsdd;
  Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
  Real tolerance=1e-6 "tolerance for RES_s (in J/kgK)";
  Integer iter = 0;
  constant Integer iter_max = 200;

algorithm
  state.phase := phase;

  if (state.phase == 2) then
    assert(T >= T_trip, "setState_Ts_error: pressure is lower than triple point pressure");
    assert(T <= T_crit, "setState_TsX_error: pressure is higher than critical pressure");
    sat := setSat_T(T=T);
    assert(s >= sat.liq.s, "setState_TsX_error: entropy is lower than saturated liquid entropy: this is single phase liquid");
    assert(s <= sat.vap.s, "setState_TsX_error: entropy is higher than saturated vapor entropy: this is single phase vapor");
  else
    if ((T <= T_crit) and (T >= T_trip)) then
      // two-phase possible, do simple check first
      sat.Tsat := T;
      tau := T_crit/sat.Tsat;
      sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
      delta := sat.liq.d/d_crit;
      f.i   := EoS.f_i(tau=tau, delta=delta);
      f.it  := EoS.f_it(tau=tau, delta=delta);
      f.r   := EoS.f_r(tau=tau, delta=delta);
      f.rt  := EoS.f_rt(tau=tau, delta=delta);
      sat.liq.s := R*(tau*(f.it + f.rt) - f.i - f.r);

      sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
      delta := sat.vap.d/d_crit;
      f.i   := EoS.f_i(tau=tau, delta=delta);
      f.it  := EoS.f_it(tau=tau, delta=delta);
      f.r   := EoS.f_r(tau=tau, delta=delta);
      f.rt  := EoS.f_rt(tau=tau, delta=delta);
      sat.vap.s := R*(tau*(f.it + f.rt) - f.i - f.r);

      if ((s > sat.liq.s - abs(0.05*sat.liq.s)) and (s < sat.vap.s + abs(0.05*sat.vap.s))) then
        // Modelica.Utilities.Streams.print("two-phase state or close to it, get saturation properties from EoS", "printlog.txt");
        sat := setSat_T(T=sat.Tsat);
      end if;

      if (s < sat.liq.s) then
        // Modelica.Utilities.Streams.print("single phase liquid: d is between dliq and rho_max", "printlog.txt");
        state.phase := 1;
        d_min := sat.liq.d;
        d_iter := sat.liq.d;
      elseif (s > sat.vap.s) then
        // Modelica.Utilities.Streams.print("single phase vapor: d is between 0 and dvap", "printlog.txt");
        state.phase := 1;
        d_max := sat.vap.d;
        d_iter := sat.vap.d/100;
      else
        // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
        state.phase := 2;
      end if;

    else
      // Modelica.Utilities.Streams.print("T>T_crit or T<T_trip, only single phase possible, do not change dmin and dmax", "printlog.txt");
      state.phase := 1;
      d_iter := d_crit/100;
    end if;
  end if;

  // phase and region determination finished !

  if (state.phase == 2) then
    // force two-phase, SaturationProperties are already known
    state.T := T;
    state.s := s;
    state.p := sat.psat;
    x := (s - sat.liq.s)/(sat.vap.s - sat.liq.s);
    state.d := 1/(1/sat.liq.d + x*(1/sat.vap.d - 1/sat.liq.d));
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
  else
    // force single-phase

    // calculate RES_s
    delta := d_iter/d_crit;
    f.i  := EoS.f_i(delta=delta, tau=tau);
    f.it := EoS.f_it(delta=delta, tau=tau);
    f.r  := EoS.f_r(delta=delta, tau=tau);
    f.rt := EoS.f_rt(delta=delta, tau=tau);
    RES_s := R*(tau*(f.it + f.rt) - f.i - f.r) - s;

    while ((abs(RES_s) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;

      // calculate gradient with respect to density
      f.rd := EoS.f_rd(delta=delta, tau=tau);
      f.rtd := EoS.f_rtd(delta=delta, tau=tau);
      dsdd := R/d_iter*(-(1+delta*f.rd-delta*tau*f.rtd));

      // print for debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("d_iter=" + String(d_iter) + " and dsdd=" + String(dsdd), "printlog.txt");

      // calculate better d_iter
      d_iter := d_iter - gamma/dsdd*RES_s;

      // check bounds
      if (d_iter<d_min) or (d_iter>d_max) then
        // Modelica.Utilities.Streams.print("d_iter out of bounds, fallback to Ridders' method, step=" + String(iter) + ", d_iter=" + String(d_iter), "printlog.txt");
        // calculate RES_s for d_min
        delta := d_min/d_crit;
        f.i  := EoS.f_i(delta=delta, tau=tau);
        f.it := EoS.f_it(delta=delta, tau=tau);
        f.r  := EoS.f_r(delta=delta, tau=tau);
        f.rt := EoS.f_rt(delta=delta, tau=tau);
        RES_min := R*(tau*(f.it + f.rt) - f.i - f.r) - s;
        // calculate RES_s for d_max
        delta := d_max/d_crit;
        f.i  := EoS.f_i(delta=delta, tau=tau);
        f.it := EoS.f_it(delta=delta, tau=tau);
        f.r  := EoS.f_r(delta=delta, tau=tau);
        f.rt := EoS.f_rt(delta=delta, tau=tau);
        RES_max := R*(tau*(f.it + f.rt) - f.i - f.r) - s;
        // calculate RES_s for d_med
        d_med := (d_max+1*d_min)/2;
        delta := d_med/d_crit;
        f.i  := EoS.f_i(delta=delta, tau=tau);
        f.it := EoS.f_it(delta=delta, tau=tau);
        f.r  := EoS.f_r(delta=delta, tau=tau);
        f.rt := EoS.f_rt(delta=delta, tau=tau);
        RES_med := R*(tau*(f.it + f.rt) - f.i - f.r) - s;
        // Ridders' method
        d_iter := d_med + (d_med-d_min)*sign(RES_min-RES_max)*RES_med/sqrt(RES_med^2-RES_min*RES_max);
        // calculate new RES_s
        delta := d_iter/d_crit;
        f.i  := EoS.f_i(delta=delta, tau=tau);
        f.it := EoS.f_it(delta=delta, tau=tau);
        f.r  := EoS.f_r(delta=delta, tau=tau);
        f.rt := EoS.f_rt(delta=delta, tau=tau);
        RES_s := R*(tau*(f.it + f.rt) - f.i - f.r) - s;
        // thighten the bounds
        if (RES_s*RES_med<=0) then
          // opposite sign, d_med and d_iter bracket the root
          d_min := min(d_med,d_iter);
          d_max := max(d_med,d_iter);
        else
          if (RES_s*RES_min<0) then
            d_max := d_iter;
          elseif (RES_s*RES_max<0) then
            d_min := d_iter;
          else
            assert(false,"never get here");
          end if;
        end if;
        // Modelica.Utilities.Streams.print("Ridders' method: new d_min=" + String(d_min) + ", new d_max=" + String(d_max), "printlog.txt");
      else
        // use d_iter from Newton
        // calculate new RES_s
        delta := d_iter/d_crit;
        f.i  := EoS.f_i(delta=delta, tau=tau);
        f.it := EoS.f_it(delta=delta, tau=tau);
        f.r  := EoS.f_r(delta=delta, tau=tau);
        f.rt := EoS.f_rt(delta=delta, tau=tau);
        RES_s := R*(tau*(f.it + f.rt) - f.i - f.r) - s;
        end if;
    end while;
    // Modelica.Utilities.Streams.print("setState_Ts total iteration steps " + String(iter) + " for T=" + String(T) + " and s=" + String(s), "printlog.txt");
    // Modelica.Utilities.Streams.print(" ", "printlog.txt");
    assert(iter<iter_max, "setState_Ts did not converge, input was T=" + String(T) + " and s=" + String(s));

    state.T := T;
    state.s :=s;
    state.d := d_iter;
    state.p := state.d*T*R*(1+delta*f.rd);
    state.h := state.T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
    state.u := state.T*R*(tau*(f.it+f.rt));
  end if;

end setState_Ts;
