within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_h_liq
  "ancillary iterative function: calculate saturation temperature for a given enthalpy of saturated liquid"

  input Modelica.Units.SI.SpecificEnthalpy h;
  output Modelica.Units.SI.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/fluidConstants[1].molarMass "specific gas constant";
  constant Real h_scale=R_s*T_crit "scale for the residual";

  // the bracket ends at 0.99*T_crit: closer to the critical point the ancillary liquid density
  // (fitted to slightly different critical values) is not accurate enough for a monotonic enthalpy
  constant Temperature T_hi=0.99*T_crit;
  constant EoS.HelmholtzDerivs fl_trip = EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T_trip), T=T_trip, order=1);
  constant EoS.HelmholtzDerivs fl_hi = EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T_hi), T=T_hi, order=1);
  constant SpecificEnthalpy hl_trip = EoS.h(fl_trip);
  constant SpecificEnthalpy hl_hi = EoS.h(fl_hi);

  Temperature T1=T_trip;
  Temperature T2=T_hi;
  Temperature T3;
  Temperature T4;

  EoS.HelmholtzDerivs f3;
  EoS.HelmholtzDerivs f4;

  SpecificEnthalpy R1 "residual of T1";
  SpecificEnthalpy R2 "residual of T2";
  SpecificEnthalpy R3 "residual of T3";
  SpecificEnthalpy R4 "residual of T4";

  constant Real tolerance=1e-6 "tolerance for the residual (relative to h_scale) and for T";
  Boolean converged=false;
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
  if (h<hl_hi) and (h>hl_trip) then
    R1 := hl_trip-h;
    R2 := hl_hi-h;
    while (iter<iter_max) and not converged loop
      iter:=iter+1;
      T3 := (T1+T2)/2;
      f3 := EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T3), T=T3, order=1);
      R3 := EoS.h(f3)-h;
      // calculate better T from Ridders' method
      T4  := T3 + (T3 - T1)*sign(R1-R2)*R3/sqrt(R3*R3 - R1*R2);
      f4 := EoS.setHelmholtzDerivs(d=Ancillary.bubbleDensity_T(T4), T=T4, order=1);
      R4 := EoS.h(f4)-h;
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
        assert(false, "Ancillary.saturationTemperature_h_liq: this should not happen. Input was h=" + String(h));
      end if;
      converged := abs(R4)<=tolerance*h_scale or abs(T1-T2)<=tolerance;
    end while;
    assert(converged, "Ancillary.saturationTemperature_h_liq did not converge, input was h=" + String(h));
    T := T4;
  elseif (h>=hl_hi) then
    T := T_crit;
  else
    T := T_trip;
  end if;

end saturationTemperature_h_liq;
