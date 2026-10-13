within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.Ancillary;
function setSat_T_start
  "start values for setSat_T: reduced saturated liquid and vapour density from the ancillary equations"
  input Temperature T;
  output Real delta_liq(unit="1") "reduced liquid density, above the solution";
  output Real delta_vap(unit="1") "reduced vapour density, below the solution";

protected
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant Density d_crit=fluidConstants[1].molarMass/fluidConstants[1].criticalMolarVolume;
  constant Density d_reducing=helmholtzCoefficients.d_reducing "reducing density of the EoS";
  constant Real delta_crit=d_crit/d_reducing "critical density, reduced with the EoS reducing density";
  constant Real delta_liq_max=1.1*fluidLimits.DMAX/d_reducing "upper bound for reduced liquid density";
  constant Real delta_vap_min=fluidLimits.DMIN/d_reducing "lower bound for reduced vapour density";
  Temperature T_anc "temperature at which the ancillary equations are evaluated";

algorithm
  // Akasaka (2008) requires the liquid guess above and the vapour guess below the solution
  if (T>T_crit*(1-1e-2)) then
    // close to the critical point, evaluate the ancillary equations in reduced variables, so that their critical point
    // coincides with the critical point of the EoS (they may be fitted to slightly different critical values)
    T_anc := T*ancillaryCoefficients.T_crit/T_crit;
    delta_liq := 1.02*Ancillary.bubbleDensity_T(T=T_anc)*d_crit/ancillaryCoefficients.d_crit/d_reducing;
    delta_vap := 0.98*Ancillary.dewDensity_T(T=T_anc)*d_crit/ancillaryCoefficients.d_crit/d_reducing;
  elseif (T>=fluidLimits.TMIN) then
    delta_liq := 1.02*Ancillary.bubbleDensity_T(T=T)/d_reducing;
    delta_vap := 0.98*Ancillary.dewDensity_T(T=T)/d_reducing;
  else
    // below TMIN the ancillary equations are extrapolated and not reliable, use safe (but less accurate) guesses
    delta_liq := delta_liq_max;
    delta_vap := 0.01*Ancillary.dewDensity_T(T=T)/d_reducing;
  end if;
  delta_liq := min(max(delta_liq, delta_crit), delta_liq_max);
  delta_vap := min(max(delta_vap, delta_vap_min), delta_crit);
end setSat_T_start;
