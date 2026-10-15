within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_derd_T_cache_der "time derivative of specificEnthalpy_derd_T_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real dhdT_der "time derivative of (dh/dd)@T";
protected
  EoS.HelmholtzDerivs f3 "Helmholtz derivatives up to third order at (d, T)";
  Real A "(dp/dd)@T";
  Real B "(dp/dT)@d";
  Real C "(dh/dd)@T";
  Real D "(dh/dT)@d";
  Real S_d "(ds/dd)@T";
  Real S_T "(ds/dT)@d";
  Real A_t;
  Real B_t;
  Real C_t;
  Real D_t;
  Real U_T_t "time derivative of (du/dT)@d";
  Real S_d_t;
  Real S_T_t;
  Real vl "specific volume of saturated liquid";
  Real vv "specific volume of saturated vapour";
  Real hl "specific enthalpy of saturated liquid";
  Real hv "specific enthalpy of saturated vapour";
  Real ul "specific inner energy of saturated liquid";
  Real uv "specific inner energy of saturated vapour";
  Real sl "specific entropy of saturated liquid";
  Real sv "specific entropy of saturated vapour";
  Real P1 "dp_sat/dT along the saturation line";
  Real P2 "d2p_sat/dT2 along the saturation line";
  Real rhol_T "d(rho_liq)/dT along saturation";
  Real rhov_T "d(rho_vap)/dT along saturation";
  Real hl_T "d(h_liq)/dT along saturation";
  Real hv_T "d(h_vap)/dT along saturation";
  Real ul_T "d(u_liq)/dT along saturation";
  Real uv_T "d(u_vap)/dT along saturation";
  Real sl_T "d(s_liq)/dT along saturation";
  Real sv_T "d(s_vap)/dT along saturation";
  Real vl_T "d(v_liq)/dT along saturation";
  Real vv_T "d(v_vap)/dT along saturation";
  Real K;
  Real K_T;
algorithm
  if phase==1 then
    f3 := EoS.setHelmholtzDerivs(d=d, T=T, order=3);
    A := EoS.dpdT(f3);
    B := EoS.dpTd(f3);
    C := EoS.dhdT(f3);
    D := EoS.dhTd(f3);
    S_d := EoS.dsdT(f3);
    S_T := EoS.dsTd(f3);
    A_t := EoS.d2pd2T(f3)*d_der + EoS.d2pTd(f3)*T_der;
    B_t := EoS.d2pTd(f3)*d_der + EoS.d2pT2d(f3)*T_der;
    C_t := EoS.d2hd2T(f3)*d_der + EoS.d2hTd(f3)*T_der;
    D_t := EoS.d2hTd(f3)*d_der + EoS.d2hT2d(f3)*T_der;
    U_T_t := EoS.d2uTd(f3)*d_der + EoS.d2uT2d(f3)*T_der;
    S_d_t := EoS.d2sd2T(f3)*d_der + EoS.d2sTd(f3)*T_der;
    S_T_t := EoS.d2sTd(f3)*d_der + EoS.d2sT2d(f3)*T_der;
    dhdT_der := C_t;
  elseif phase==2 then
    vl := 1/cache.sat.liq.d;
    vv := 1/cache.sat.vap.d;
    hl := cache.sat.liq.h;
    hv := cache.sat.vap.h;
    ul := cache.sat.liq.u;
    uv := cache.sat.vap.u;
    sl := cache.sat.liq.s;
    sv := cache.sat.vap.s;
    P1 := (sv - sl)/(vv - vl) "Clausius-Clapeyron";
    rhol_T := (P1 - EoS.dpTd(cache.fl))/EoS.dpdT(cache.fl);
    rhov_T := (P1 - EoS.dpTd(cache.fv))/EoS.dpdT(cache.fv);
    hl_T := EoS.dhTd(cache.fl) + EoS.dhdT(cache.fl)*rhol_T;
    hv_T := EoS.dhTd(cache.fv) + EoS.dhdT(cache.fv)*rhov_T;
    ul_T := EoS.duTd(cache.fl) + EoS.dudT(cache.fl)*rhol_T;
    uv_T := EoS.duTd(cache.fv) + EoS.dudT(cache.fv)*rhov_T;
    sl_T := EoS.dsTd(cache.fl) + EoS.dsdT(cache.fl)*rhol_T;
    sv_T := EoS.dsTd(cache.fv) + EoS.dsdT(cache.fv)*rhov_T;
    vl_T := -rhol_T*vl^2;
    vv_T := -rhov_T*vv^2;
    P2 := ((sv_T - sl_T)*(vv - vl) - (sv - sl)*(vv_T - vl_T))/(vv - vl)^2;
    // two phase: dhdT = -K/d^2, K = (h'-h")/(v'-v") depends on T only
    K := (hl - hv)/(vl - vv);
    K_T := ((hl_T - hv_T)*(vl - vv) - (hl - hv)*(vl_T - vv_T))/(vl - vv)^2;
    dhdT_der := 2*K/d^3*d_der - K_T/d^2*T_der;
  else
    assert(false, "specificEnthalpy_derd_T_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEnthalpy_derd_T_cache\">specificEnthalpy_derd_T_cache</a>,
used for analytic Jacobians. In the two-phase region (dh/dd)@T = -1/d^2*(h'-h\")/(v'-v\"), its derivative needs the first derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end specificEnthalpy_derd_T_cache_der;
