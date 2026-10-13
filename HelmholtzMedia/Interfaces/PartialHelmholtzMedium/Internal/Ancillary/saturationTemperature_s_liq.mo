within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function saturationTemperature_s_liq
  "ancillary iterative function: calculate saturation temperature for a given entropy of saturated liquid"

  input Modelica.Units.SI.SpecificEntropy s;
  output Modelica.Units.SI.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/fluidConstants[1].molarMass "specific gas constant";
  constant Real s_scale=R_s "scale for the residual";

  // the bracket ends at 0.99*T_crit: closer to the critical point the ancillary liquid density
  // (fitted to slightly different critical values) is not accurate enough for a monotonic entropy
  constant Temperature T_hi=0.99*T_crit;
  constant EoS.HelmholtzDerivs fl_trip = EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T_trip), T=T_trip, order=1);
  constant EoS.HelmholtzDerivs fl_hi = EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T_hi), T=T_hi, order=1);
  constant SpecificEntropy sl_trip = EoS.s(fl_trip);
  constant SpecificEntropy sl_hi = EoS.s(fl_hi);

  Temperature T1=T_trip;
  Temperature T2=T_hi;
  Temperature T3;
  Temperature T4;

  EoS.HelmholtzDerivs f3;
  EoS.HelmholtzDerivs f4;

  SpecificEntropy R1 "residual of T1";
  SpecificEntropy R2 "residual of T2";
  SpecificEntropy R3 "residual of T3";
  SpecificEntropy R4 "residual of T4";

  constant Real tolerance=1e-6 "tolerance for the residual (relative to s_scale) and for T";
  Boolean converged=false;
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
  if (s<sl_hi) and (s>sl_trip) then
    R1 := sl_trip-s;
    R2 := sl_hi-s;
    while (iter<iter_max) and not converged loop
      iter:=iter+1;
      T3 := (T1+T2)/2;
      f3 := EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T3), T=T3, order=1);
      R3 := EoS.s(f3)-s;
      // calculate better T from Ridders' method
      T4  := T3 + (T3 - T1)*sign(R1-R2)*R3/sqrt(R3*R3 - R1*R2);
      f4 := EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T4), T=T4, order=1);
      R4 := EoS.s(f4)-s;
      if (R4*R3<=0) then
        T1 := T3;
        R1 := R3;
        T2 := T4;
        R2 := R4;
      elseif (R4*R1<0) then
        T2 := T4;
        R2 := R4;
      elseif (R4*R2<0) then
        T1 := T4;
        R1 := R4;
      else
        assert(false, "Ancillary.saturationTemperature_s_liq: this should not happen. Input was s=" + String(s));
      end if;
      converged := abs(R4)<=tolerance*s_scale or abs(T1-T2)<=tolerance;
    end while;
    assert(converged, "Ancillary.saturationTemperature_s_liq did not converge, input was s=" + String(s));
    T := T4;
  elseif (s>=sl_hi) then
    T := T_crit;
  else
    T := T_trip;
  end if;

end saturationTemperature_s_liq;
