within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setState_pd "Return thermodynamic state as function of (p, d)"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input Density d "Density";
  input FixedPhase phase=0 "2 for two-phase, 1 for one-phase, 0 if not known";
  output ThermodynamicState state "thermodynamic state record";

protected
  MolarMass MM = fluidConstants[1].molarMass;
  SpecificHeatCapacity R=Modelica.Constants.R/MM "specific gas constant";
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real delta(unit="1")=d/d_crit "reduced density";
  Real tau(unit="1") "inverse reduced temperature";
  HelmholtzDerivs f;

  AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  SaturationProperties sat;
  MassFraction x "vapour quality";
  Temperature Tmin=fluidLimits.TMIN;
  Temperature Tmax=fluidLimits.TMAX;
  Real tolerance=1e-9 "relative Tolerance for Density";

algorithm
  state.phase := phase;

  if (state.phase == 2) then
    assert(p >= p_trip, "setState_pdX_error: pressure is lower than triple point pressure");
    assert(p <= p_crit, "setState_pdX_error: pressure is higher than critical pressure");
    sat := setSat_p(p=p);
    assert(d <= sat.liq.d, "setState_pdX_error: density is higher than saturated liquid density: this is single phase liquid");
    assert(d >= sat.vap.d, "setState_pdX_error: density is lower than saturated vapor density: this is single phase vapor");
  else
    if ((p <= p_crit) and (p >= p_trip)) then
      // two-phase possible, do simple check first
      sat.Tsat := saturationTemperature(p=p);
      sat.liq.d := bubbleDensity_T_ANC(T=sat.Tsat);
      sat.vap.d := dewDensity_T_ANC(T=sat.Tsat);
      // Modelica.Utilities.Streams.print("setState_pd: sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d) + ", simple check only");

      if ((d < sat.liq.d + abs(0.05*sat.liq.d)) and (d > sat.vap.d - abs(0.05*sat.vap.d))) then
        // Modelica.Utilities.Streams.print("setState_pd: p = " + String(p) + "d = " + String(d) + ", two-phase state or close to it");
        // get saturation properties from EoS, use Tsat as starting value
        sat := setSat_p(p=p, T_guess=sat.Tsat);
        // Modelica.Utilities.Streams.print("setState_pd: sat.liq.d=" + String(sat.liq.d) + " sat.vap.d=" + String(sat.vap.d) + ", from EoS");
      end if;

      if (d > sat.liq.d) then
        state.phase := 1; // single phase liquid
        Tmax := sat.Tsat;
        // Tmin := Tsat(d);
      elseif (d < sat.vap.d) then
        state.phase := 1; // single phase vapor
        Tmin := sat.Tsat;
      else
        state.phase := 2; // two-phase, all properties can be calculated from sat record
      end if;

    else
      // p>p_crit or p<p_trip, only single phase possible, do not change Tmin and Tmax
      state.phase := 1;
    end if;
  end if;

  state.p := p;
  state.d := d;
  if (state.phase == 2) then
    // force two-phase, SaturationProperties are already known
    x := (1/d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
    state.T := sat.Tsat;
    state.h := sat.liq.h + x*(sat.vap.h - sat.liq.h);
    state.u := sat.liq.u + x*(sat.vap.u - sat.liq.u);
    state.s := sat.liq.s + x*(sat.vap.s - sat.liq.s);
  else
    // force single-phase
    state.T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
          function HelmholtzMedia.Interfaces.PartialHelmholtzMedium.setState_pd_RES(
            p=p,
            d=d,
            phase=1),
          u_min=Tmin,
          u_max=Tmax,
          tolerance=tolerance);

    tau := T_crit/state.T;
    delta := state.d/d_crit;

    f.i   := f_i(tau=tau, delta=delta);
    f.it  := f_it(tau=tau, delta=delta);
    f.r   := f_r(tau=tau, delta=delta);
    f.rt  := f_rt(tau=tau, delta=delta);
    f.rd  := f_rd(tau=tau, delta=delta);
    state.h := state.T*R*(tau*(f.it + f.rt) + (1+delta*f.rd));
    state.u := state.T*R*(tau*(f.it+f.rt));
    state.s :=         R*(tau*(f.it+f.rt) - (f.i+f.r));
  end if;

end setState_pd;
