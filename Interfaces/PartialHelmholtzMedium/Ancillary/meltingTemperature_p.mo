within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function meltingTemperature_p
  "ancillary iterative function: calculate melting temperature for a given pressure by iterating the ancillary function"
  input AbsolutePressure p;
  output Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_max=fluidLimits.TMAX;

  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_max=fluidLimits.PMAX;

  Temperature T1=0.98*T_trip "low temperature";
  Temperature T2=T_max "high temperature";
  Temperature T3 "intermediate temperature";
  Temperature T4 "new temperature";

  Real R1 "residual of T1";
  Real R2 "residual of T2";
  Real R3 "residual of T3";
  Real R4= Modelica.Constants.inf "residual of T4";

  constant Real tolerance=1e-6 "relative tolerance for RES_p";
  Integer iter=0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print("Ancillary.meltingTemperature_p, p=" + String(p), "printlog.txt");

  if (p<p_max) and (p>p_trip) then
    R1 := Ancillary.meltingPressure_T(T1)-p;
    R2 := Ancillary.meltingPressure_T(T2)-p;
    if (R1*R2<0) then
      while (abs(R4/p)>tolerance) and (iter<iter_max) and (abs(T1-T2)>tolerance) loop
        iter:=iter+1;
        T3 := (T1+T2)/2;
        R3 := Ancillary.meltingPressure_T(T3)-p;
        // calculate better T from Ridder's method
        T4  := T3 + (T3 - T1)*sign(R1-R2)*R3/sqrt(R3*R3 - R1*R2);
        R4 := Ancillary.meltingPressure_T(T4)-p;
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
            assert(false, "Ancillary.meltingTemperature_p: this should not happen");
          end if;
        end if;
        // Modelica.Utilities.Streams.print("Ridders' method: new brackets T1=" + String(T1) + " and T2=" + String(T2), "printlog.txt");
      end while;
      assert(iter<iter_max, "meltingTemperature_p did not converge, input was p=" + String(p));
      // Modelica.Utilities.Streams.print("Ancillary.meltingTemperature_p total iteration steps " + String(iter) + " for p=" + String(p), "printlog.txt");
      // Modelica.Utilities.Streams.print(" ", "printlog.txt");
      T := T4;
    else
      if (abs(R1/p)<tolerance) then
        T:= T1;
      elseif (abs(R2/p)<tolerance) then
        T:=T2;
      else
        assert(false, "Ancillary.meltingTemperature_p: T1=" + String(T1) + " and T2=" + String(T2) + " did not bracket the root");
      end if;
    end if;

  elseif (p<=p_trip) then
    T := T_trip;
  else
    T := T_max;
  end if;

  annotation (inverse(p=meltingPressure_T(T=T)));
end meltingTemperature_p;
