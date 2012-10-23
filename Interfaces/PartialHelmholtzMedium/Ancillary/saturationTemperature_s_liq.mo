within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function saturationTemperature_s_liq
  "ancillary iterative function: calculate saturation temperature for a given entropy by iterating the ancillary function"

  input Modelica.SIunits.SpecificEntropy s;
  output Modelica.SIunits.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant SpecificEntropy s_crit=fluidConstants[1].SCRIT0;

  constant Density dl_trip = Ancillary.bubbleDensity_T(T_trip);
  constant EoS.HelmholtzDerivs fl_trip = EoS.setHelmholtzDerivsFirst(d=dl_trip, T=T_trip);
  constant SpecificEntropy sl_trip = EoS.s(fl_trip);

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

  SpecificEntropy R1 "residual of T1";
  SpecificEntropy R2 "residual of T2";
  SpecificEntropy R3 "residual of T3";
  SpecificEntropy R4= Modelica.Constants.inf "residual of T4";

  constant Real tolerance=1e-6 "relative tolerance for RES";
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print("Ancillary.saturationTemperature_h, h=" + String(h), "printlog.txt");
  // assert(s<=s_crit+tolerance,"bubble entropy cannot be higher than critical entropy, invalid input");

  if (s<0.98*s_crit) and (s>sl_trip) then
    // Modelica.Utilities.Streams.print("s=" + String(s,significantDigits=9) + "<s_crit: liquid side", "printlog.txt");
    R1 := sl_trip-s;
    R2 := s_crit-s;
    if (R1*R2<0) then
      while (abs(R4/s)>tolerance) and (iter<iter_max) and (abs(T1-T2)>tolerance) loop
        iter:=iter+1;
        T3 := (T1+T2)/2;
        d3 := Ancillary.bubbleDensity_T(T3);
        f3 := EoS.setHelmholtzDerivsFirst(d=d3,T=T3);
        R3 := EoS.s(f3)-s;
        // caclutate better T from Ridder's method
        T4  := T3 + (T3 - T1)*sign(R1-R2)*R3/sqrt(R3*R3 - R1*R2);
        d4 := Ancillary.bubbleDensity_T(T4);
        f4 := EoS.setHelmholtzDerivsFirst(d=d4,T=T4);
        R4 := EoS.s(f4)-s;
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
            assert(false, "Ancillary.saturationTemperature_s (liquid side): this should not happen");
          end if;
        end if;
        // Modelica.Utilities.Streams.print("Ridders' method: new brackets T1=" + String(T1) + " and T2=" + String(T2), "printlog.txt");
      end while;
      assert(iter<iter_max, "saturationTemperature_s_liq did not converge, input was s_liq=" + String(s));
      // Modelica.Utilities.Streams.print("Ancillary.saturationTemperature_s_liq total iteration steps " + String(iter) + " for s_liq=" + String(s), "printlog.txt");
      // Modelica.Utilities.Streams.print(" ", "printlog.txt");
      T := T4;
    else
      if (abs(R1/s)<tolerance) then
        T:= T1;
      elseif (abs(R2/s)<tolerance) then
        T:=T2;
      else
        assert(false, "Ancillary.saturationTemperature_s (liquid side): T1=" + String(T1) + " and T2=" + String(T2) + " did not bracket the root");
      end if;
    end if;

  elseif (s>=0.98*s_crit) then
    T := T_crit;
  elseif (s<=sl_trip) then
    T := T_trip;
  else
    assert(false, "Ancillary.saturationTemperature_s (liquid side): this should also not happen");
  end if;

end saturationTemperature_s_liq;
