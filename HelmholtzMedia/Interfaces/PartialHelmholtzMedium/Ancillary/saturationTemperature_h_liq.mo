within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_h_liq
  "ancillary iterative function: calculate saturation temperature for a given enthalpy by iterating the ancillary function"

  input Modelica.Units.SI.SpecificEnthalpy h;
  output Modelica.Units.SI.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant SpecificEnthalpy h_crit=fluidConstants[1].HCRIT0;

  constant Density dl_trip = Ancillary.bubbleDensity_T(T_trip);
  constant EoS.HelmholtzDerivs fl_trip = EoS.setHelmholtzDerivsFirst(d=dl_trip, T=T_trip);
  constant SpecificEnthalpy hl_trip = EoS.h(fl_trip);

  Temperature T1=0.98*T_trip;
  Temperature T2=T_crit;
  Temperature T3;
  Temperature T4;

  Density d1;
  Density d2;
  Density d3;
  Density d4;

  EoS.HelmholtzDerivs f1;
  EoS.HelmholtzDerivs f2;
  EoS.HelmholtzDerivs f3;
  EoS.HelmholtzDerivs f4;

  SpecificEnthalpy R1 "residual of T1";
  SpecificEnthalpy R2 "residual of T2";
  SpecificEnthalpy R3 "residual of T3";
  SpecificEnthalpy R4= Modelica.Constants.inf "residual of T4";

  constant Real tolerance=1e-6 "relative tolerance for RES";
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print("Ancillary.saturationTemperature_h, h=" + String(h), "printlog.txt");
  // assert(h<=h_crit+tolerance,"bubble enthalpy cannot be higher than critical enthalpy, invalid input");

  if (h<0.98*h_crit) and (h>hl_trip) then
    // Modelica.Utilities.Streams.print("h<h_crit: liquid side", "printlog.txt");
    R1 := hl_trip-h;
    R2 := h_crit-h;
    if (R1*R2<0) then
      while (abs(R4/h)>tolerance) and (iter<iter_max) and (abs(T1-T2)>tolerance) loop
        iter:=iter+1;
        T3 := (T1+T2)/2;
        d3 := Ancillary.bubbleDensity_T(T3);
        f3 := EoS.setHelmholtzDerivsFirst(d=d3,T=T3);
        R3 := EoS.h(f3)-h;
        // caclutate better T from Ridder's method
        T4  := T3 + (T3 - T1)*sign(R1-R2)*R3/sqrt(R3*R3 - R1*R2);
        d4 := Ancillary.bubbleDensity_T(T4);
        f4 := EoS.setHelmholtzDerivsFirst(d=d4,T=T4);
        R4 := EoS.h(f4)-h;
        // Modelica.Utilities.Streams.print("Ridders' method: current residuals: R1=" + String(R1) + ", R2=" + String(R2) + ", R3=" + String(R3) + ", R4=" + String(R4), "printlog.txt");
        if (R4*R3<=0) then
          T1 := T3;
          R1 := R3;
          T2 := T4;
          R2 := R4;
        else
          if (R4*R1<0) then
            T2 := T4;
            R2 := R4;
          elseif (R4*R2<0) then
            T1 := T4;
            R1 := R4;
          else
            assert(false, "Ancillary.saturationTemperature_h (liquid side): this should not happen. Input was h=" + String(h));
          end if;
        end if;
        // Modelica.Utilities.Streams.print("Ridders' method: new brackets T1=" + String(T1) + " and T2=" + String(T2), "printlog.txt");
      end while;
      assert(iter<iter_max, "saturationTemperature_h_liq did not converge, input was h_liq=" + String(h));
      // Modelica.Utilities.Streams.print("Ancillary.saturationTemperature_h_liq total iteration steps " + String(iter) + " for h_liq=" + String(h), "printlog.txt");
      // Modelica.Utilities.Streams.print(" ", "printlog.txt");
      T := T4;
    else
      if (abs(R1/h)<tolerance) then
        T:= T1;
      elseif (abs(R2/h)<tolerance) then
        T:=T2;
      else
        assert(false, "Ancillary.saturationTemperature_h (liquid side): T1=" + String(T1) + " and T2=" + String(T2) + " did not bracket the root");
      end if;
    end if;

  elseif (h>=0.98*h_crit) then
    T := T_crit;
  elseif (h<=hl_trip) then
    T := T_trip;
  else
    assert(false, "Ancillary.saturationTemperature_h (liquid side): this should also not happen. Input was h=" + String(h));
  end if;

end saturationTemperature_h_liq;
