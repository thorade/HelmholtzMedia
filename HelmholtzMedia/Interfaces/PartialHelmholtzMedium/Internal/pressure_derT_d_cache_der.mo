within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function pressure_derT_d_cache_der "time derivative of pressure_derT_d_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real dpTd_der "time derivative of (dp/dT)@d";
protected
  EoS.HelmholtzDerivs f3 "Helmholtz derivatives up to third order at (d, T)";
  Real P1 "dp_sat/dT along the saturation line (Clausius-Clapeyron)";
  Real vl "specific volume of saturated liquid";
  Real vv "specific volume of saturated vapour";
  Real hl "specific enthalpy of saturated liquid";
  Real hv "specific enthalpy of saturated vapour";
  Real rhol_T "d(rho_liq)/dT along saturation";
  Real rhov_T "d(rho_vap)/dT along saturation";
  Real hl_T "d(h_liq)/dT along saturation";
  Real hv_T "d(h_vap)/dT along saturation";
  Real vl_T "d(v_liq)/dT along saturation";
  Real vv_T "d(v_vap)/dT along saturation";
  Real sl_T "d(s_liq)/dT along saturation";
  Real sv_T "d(s_vap)/dT along saturation";
  Real P2 "d2p_sat/dT2 along the saturation line";
algorithm
  if phase==1 then
    f3 := EoS.setHelmholtzDerivs(d=d, T=T, order=3);
    dpTd_der := EoS.d2pTd(f3)*d_der + EoS.d2pT2d(f3)*T_der;
  elseif phase==2 then
    vl := 1/cache.sat.liq.d;
    vv := 1/cache.sat.vap.d;
    hl := cache.sat.liq.h;
    hv := cache.sat.vap.h;
    P1 := (cache.sat.vap.s - cache.sat.liq.s)/(vv - vl);
    rhol_T := (P1 - EoS.dpTd(cache.fl))/EoS.dpdT(cache.fl);
    rhov_T := (P1 - EoS.dpTd(cache.fv))/EoS.dpdT(cache.fv);
    hl_T := EoS.dhTd(cache.fl) + EoS.dhdT(cache.fl)*rhol_T;
    hv_T := EoS.dhTd(cache.fv) + EoS.dhdT(cache.fv)*rhov_T;
    vl_T := -rhol_T*vl^2;
    vv_T := -rhov_T*vv^2;
    sl_T := (hl_T - vl*P1)/T "ds = (dh - v*dp)/T along saturation";
    sv_T := (hv_T - vv*P1)/T;
    P2 := ((sv_T - sl_T)*(vv - vl) - (cache.sat.vap.s - cache.sat.liq.s)*(vv_T - vl_T))/(vv - vl)^2;
    // two phase: (dp/dT)@d = dp_sat/dT, a function of T only
    dpTd_der := P2*T_der;
  else
    assert(false, "pressure_derT_d_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.pressure_derT_d_cache\">pressure_derT_d_cache</a>,
used for analytic Jacobians. In the two-phase region (dp/dT)@d is the slope of the vapour pressure curve; its derivative needs the first derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end pressure_derT_d_cache_der;
