within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_pd "Return thermodynamic state as function of (p, d)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output ThermodynamicState state "thermodynamic state record";

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1") "inverse reduced temperature";
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
  constant Density dv_trip = Ancillary.dewDensity_T(T_trip);
  constant Density dl_trip = Ancillary.bubbleDensity_T(T_trip);

  EoS.HelmholtzDerivs f;
  SaturationProperties sat;
  MassFraction x "vapour quality";

  Temperature T_min=fluidLimits.TMIN;
  Temperature T_max=fluidLimits.TMAX;
  Temperature T_iter=T_crit;
  Real RES_p;
  Real dpdT;
  Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
  Real tolerance=1e-9 "relative tolerance for RES_p";
  Integer iter = 0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print(" ", "printlog.txt");
  // Modelica.Utilities.Streams.print("setState_pdX: p=" + String(p) + " and d=" + String(d), "printlog.txt");
  state.phase := phase;

  if (state.phase == 2) then
    assert(p <= p_crit, "setState_pdX_error: pressure is higher than critical pressure");
    sat := setSat_p(p=p);
    assert(d <= sat.liq.d, "setState_pdX_error: density is higher than saturated liquid density: this is single phase liquid");
    assert(d >= sat.vap.d, "setState_pdX_error: density is lower than saturated vapor density: this is single phase vapor");
  else
    if (p < p_crit) then
      // two-phase possible, do simple check first
      sat.Tsat := Ancillary.saturationTemperature_p(p=p);
      sat.liq.d := Ancillary.bubbleDensity_T(T=sat.Tsat);
      sat.vap.d := Ancillary.dewDensity_T(T=sat.Tsat);
      // Modelica.Utilities.Streams.print("setState_pd: sat.Tsat=" + String(sat.Tsat) + " and sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d) + ", simple check only", "printlog.txt");

      if ((d < sat.liq.d + abs(0.05*sat.liq.d)) and (d > sat.vap.d - abs(0.05*sat.vap.d))) or (p<300*p_trip) or (p>0.98*p_crit) then
        // Modelica.Utilities.Streams.print("setState_pd: p = " + String(p) + "d = " + String(d) + ", two-phase state or close to it", "printlog.txt");
        // get saturation properties from EoS
        sat := setSat_p(p=p);
      end if;

      // Modelica.Utilities.Streams.print("setState_pd: phase boundary determined: sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d), "printlog.txt");
      if (d > sat.liq.d) then
        // Modelica.Utilities.Streams.print("single phase liquid", "printlog.txt");
        state.phase := 1;
        T_max := sat.Tsat;
        T_iter:= Ancillary.saturationTemperature_d(d=d); // look at subcritical isobars in T,d-Diagram !!
        T_min := 0.98*T_iter;
      elseif (d < sat.vap.d) then
        // Modelica.Utilities.Streams.print("single phase vapour", "printlog.txt");
        state.phase := 1;
        T_min := sat.Tsat;
        T_iter:= sat.Tsat;
      else
        // Modelica.Utilities.Streams.print("two-phase, all properties can be calculated from sat record", "printlog.txt");
        state.phase := 2;
      end if;

    elseif (p >= p_crit) then
      // Modelica.Utilities.Streams.print("p>p_crit , only single phase possible", "printlog.txt");
      state.phase := 1;
    else
      assert(false, "setState_pd: this should not happen, check p");
    end if;
  end if;

  // phase and region determination finished !

  if (state.phase == 2) then
    // force two-phase, SaturationProperties are already known
    state.p := p;
    state.d := d;
    x := (1/d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
    state.T := sat.Tsat;
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
  else
    // force single-phase

    // calculate RES_p
    tau := T_crit/T_iter;
    f.rd  := EoS.f_rd(delta=delta, tau=tau);
    RES_p := d*T_iter*R*(1+delta*f.rd) - p;

    while ((abs(RES_p/p) > tolerance) and (iter<iter_max)) loop
      iter := iter+1;

      // calculate gradients with respect to temperature
      f.rtd := EoS.f_rtd(delta=delta, tau=tau);
      dpdT := d*R*(1+delta*f.rd-delta*tau*f.rtd);

      // print for debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("T_iter=" + String(T_iter) + " and dpdT=" + String(dpdT), "printlog.txt");

      // calculate better d_iter and T_iter
      T_iter := T_iter - gamma/dpdT*RES_p;

      // check bounds
      T_iter := max(T_iter,0.98*T_min);
      T_iter := min(T_iter,1.02*T_max);

      // calculate new RES_p
      tau := T_crit/T_iter;
      f.rd  := EoS.f_rd(delta=delta, tau=tau);
      RES_p := d*T_iter*R*(1+delta*f.rd) - p;
    end while;
    // Modelica.Utilities.Streams.print("setState_pdX total iteration steps " + String(iter), "printlog.txt");
    assert(iter<iter_max, "setState_pdX did not converge, input was p=" + String(p) + " and d=" + String(d));

    state.p := p;
    state.d := d;
    state.T := T_iter;
    f.i   := EoS.f_i(tau=tau, delta=delta);
    f.it  := EoS.f_it(tau=tau, delta=delta);
    f.r   := EoS.f_r(tau=tau, delta=delta);
    f.rt  := EoS.f_rt(tau=tau, delta=delta);
    state.h := state.T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
    state.u := state.T*R*(tau*(f.it+f.rt));
    state.s :=         R*(tau*(f.it+f.rt) - (f.i+f.r));
  end if;

end setState_pd;
