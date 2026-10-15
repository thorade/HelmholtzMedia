within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificEnthalpy_derT_d_cache_der "time derivative of specificEnthalpy_derT_d_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real dhTd_der "time derivative of (dh/dT)@d";
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
  EoS.HelmholtzDerivs fl3 "Helmholtz derivatives up to third order at the saturated liquid";
  EoS.HelmholtzDerivs fv3 "Helmholtz derivatives up to third order at the saturated vapour";
  Real rhol_TT "d2(rho_liq)/dT2 along saturation";
  Real rhov_TT "d2(rho_vap)/dT2 along saturation";
  Real hl_TT "d2(h_liq)/dT2 along saturation";
  Real hv_TT "d2(h_vap)/dT2 along saturation";
  Real vl_TT "d2(v_liq)/dT2 along saturation";
  Real vv_TT "d2(v_vap)/dT2 along saturation";
  Real x "vapour quality";
  Real x_d "(dx/dd)@T";
  Real x_T "(dx/dT)@d";
  Real H_x "(dH/dx)@T";
  Real H_T "(dH/dT)@x";
algorithm
  if phase==1 then
    f3 := EoS.setHelmholtzDerivs(d=d, T=T, order=3);
    dhTd_der := EoS.d2hTd(f3)*d_der + EoS.d2hT2d(f3)*T_der;
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
    fl3 := EoS.setHelmholtzDerivs(d=cache.sat.liq.d, T=T, order=3);
    fv3 := EoS.setHelmholtzDerivs(d=cache.sat.vap.d, T=T, order=3);
    rhol_TT := (P2 - EoS.d2pT2d(fl3) - EoS.d2pTd(fl3)*rhol_T)/EoS.dpdT(fl3)
             - (P1 - EoS.dpTd(fl3))*(EoS.d2pTd(fl3) + EoS.d2pd2T(fl3)*rhol_T)/EoS.dpdT(fl3)^2;
    rhov_TT := (P2 - EoS.d2pT2d(fv3) - EoS.d2pTd(fv3)*rhov_T)/EoS.dpdT(fv3)
             - (P1 - EoS.dpTd(fv3))*(EoS.d2pTd(fv3) + EoS.d2pd2T(fv3)*rhov_T)/EoS.dpdT(fv3)^2;
    hl_TT := EoS.d2hT2d(fl3) + 2*EoS.d2hTd(fl3)*rhol_T + EoS.d2hd2T(fl3)*rhol_T^2 + EoS.dhdT(fl3)*rhol_TT;
    hv_TT := EoS.d2hT2d(fv3) + 2*EoS.d2hTd(fv3)*rhov_T + EoS.d2hd2T(fv3)*rhov_T^2 + EoS.dhdT(fv3)*rhov_TT;
    vl_TT := -rhol_TT*vl^2 + 2*rhol_T^2*vl^3;
    vv_TT := -rhov_TT*vv^2 + 2*rhov_T^2*vv^3;
    // two phase: H = dh'/dT + (x*dv"/dT + (1-x)*dv'/dT)*(h"-h')/(v'-v") + x*(dh"/dT - dh'/dT), x = (1/d - v')/(v"-v')
    // (the same expression as in specificEnthalpy_derT_d_cache, with derivatives along saturation by T)
    x := (1/d - vl)/(vv - vl);
    x_d := -1/(d^2*(vv - vl));
    x_T := (-vl_T*(vv - vl) - (1/d - vl)*(vv_T - vl_T))/(vv - vl)^2;
    H_x := (vv_T - vl_T)*(hv - hl)/(vl - vv) + (hv_T - hl_T);
    H_T := hl_TT + (x*vv_TT + (1 - x)*vl_TT)*(hv - hl)/(vl - vv)
         + (x*vv_T + (1 - x)*vl_T)*((hv_T - hl_T)*(vl - vv) - (hv - hl)*(vl_T - vv_T))/(vl - vv)^2
         + x*(hv_TT - hl_TT);
    dhTd_der := H_x*x_d*d_der + (H_T + H_x*x_T)*T_der;
  else
    assert(false, "specificEnthalpy_derT_d_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificEnthalpy_derT_d_cache\">specificEnthalpy_derT_d_cache</a>,
used for analytic Jacobians. In the two-phase region (dh/dT)@d depends on d (vapour quality) and on the saturation properties and their first derivatives along the saturation line; its time derivative needs the second derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end specificEnthalpy_derT_d_cache_der;
