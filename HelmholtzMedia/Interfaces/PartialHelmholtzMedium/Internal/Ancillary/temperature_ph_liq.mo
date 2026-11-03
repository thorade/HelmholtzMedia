within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function temperature_ph_liq
  "ancillary function: approximate temperature of (compressed) liquid for given p and h, used as start value"

  input Modelica.Units.SI.AbsolutePressure p;
  input Modelica.Units.SI.SpecificEnthalpy h;
  output Modelica.Units.SI.Temperature T;

protected
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  // same upper end as saturationTemperature_h_liq: closer to T_crit the ancillary liquid density is not accurate enough
  constant Temperature T_hi=0.99*T_crit;
  constant Integer n=16;
  // nodes dense towards T_hi, where h of the saturated liquid rises steeply
  constant Temperature[n] T_tab={T_trip + (T_hi - T_trip)*sin(0.5*Modelica.Constants.pi*(i - 1)/(n - 1)) for i in 1:n};
  constant SpecificEnthalpy[n] h_tab={EoS.h(EoS.setHelmholtzDerivs(d=bubbleDensity_T(T_tab[i]), T=T_tab[i], order=1)) for i in 1:n};

  SpecificEnthalpy hc "enthalpy of the saturated liquid at the same temperature (estimate)";
  Integer k;

algorithm
  hc := h;
  T := T_trip;
  for pass in 1:2 loop
    if (hc <= h_tab[1]) then
      T := T_trip;
    elseif (hc >= h_tab[n]) then
      T := T_hi;
    else
      k := 1;
      while (hc > h_tab[k + 1]) loop
        k := k + 1;
      end while;
      T := T_tab[k] + (T_tab[k + 1] - T_tab[k])*(hc - h_tab[k])/(h_tab[k + 1] - h_tab[k]);
    end if;
    if (pass == 1) then
      // compressed liquid: h(p,T) is approximately h_liq(T) + (p - p_sat(T))/d_liq(T)
      hc := h - max(p - saturationPressure_T(T=T), 0)/bubbleDensity_T(T=T);
    end if;
  end for;

  annotation (Documentation(info="<html>
<p>
Start value for the temperature of liquid and liquid-like supercritical states in <code>setState_phX</code>.
The enthalpy of the saturated liquid is tabulated at 16 temperatures between the triple point and 0.99&nbsp;T<sub>crit</sub>
(constants, evaluated once from the EoS and the ancillary liquid density) and inverted by linear interpolation.
In a second pass, the enthalpy is corrected for compression, h - (p - p<sub>sat</sub>(T))/d<sub>liq</sub>(T).
The function needs no evaluation of the EoS at run time.
</p>
</html>"));
end temperature_ph_liq;
