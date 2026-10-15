within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function velocityOfSound_cache_der "time derivative of velocityOfSound_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real a_der "time derivative of the velocity of sound";
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
  EoS.HelmholtzDerivs fl3 "Helmholtz derivatives up to third order at the saturated liquid";
  EoS.HelmholtzDerivs fv3 "Helmholtz derivatives up to third order at the saturated vapour";
  Real rhol_TT;
  Real rhov_TT;
  Real hl_TT;
  Real hv_TT;
  Real ul_TT;
  Real uv_TT;
  Real sl_TT;
  Real sv_TT;
  Real vl_TT;
  Real vv_TT;
  Real x "vapour quality";
  Real x_d "(dx/dd)@T";
  Real x_T "(dx/dT)@d";
  Real W;
  Real W_t;
  Real G;
  Real G_x;
  Real G_T;
  Real ddps;
  Real ddps_t;
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
    // a = sqrt(W), W = A - B*S_d/S_T
    W := A - B*S_d/S_T;
    W_t := A_t - (B_t*S_d + B*S_d_t)/S_T + B*S_d*S_T_t/S_T^2;
    a_der := W_t/(2*sqrt(W));
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
    fl3 := EoS.setHelmholtzDerivs(d=cache.sat.liq.d, T=T, order=3);
    fv3 := EoS.setHelmholtzDerivs(d=cache.sat.vap.d, T=T, order=3);
    rhol_TT := (P2 - EoS.d2pT2d(fl3) - EoS.d2pTd(fl3)*rhol_T)/EoS.dpdT(fl3)
             - (P1 - EoS.dpTd(fl3))*(EoS.d2pTd(fl3) + EoS.d2pd2T(fl3)*rhol_T)/EoS.dpdT(fl3)^2;
    rhov_TT := (P2 - EoS.d2pT2d(fv3) - EoS.d2pTd(fv3)*rhov_T)/EoS.dpdT(fv3)
             - (P1 - EoS.dpTd(fv3))*(EoS.d2pTd(fv3) + EoS.d2pd2T(fv3)*rhov_T)/EoS.dpdT(fv3)^2;
    hl_TT := EoS.d2hT2d(fl3) + 2*EoS.d2hTd(fl3)*rhol_T + EoS.d2hd2T(fl3)*rhol_T^2 + EoS.dhdT(fl3)*rhol_TT;
    hv_TT := EoS.d2hT2d(fv3) + 2*EoS.d2hTd(fv3)*rhov_T + EoS.d2hd2T(fv3)*rhov_T^2 + EoS.dhdT(fv3)*rhov_TT;
    ul_TT := EoS.d2uT2d(fl3) + 2*EoS.d2uTd(fl3)*rhol_T + EoS.d2ud2T(fl3)*rhol_T^2 + EoS.dudT(fl3)*rhol_TT;
    uv_TT := EoS.d2uT2d(fv3) + 2*EoS.d2uTd(fv3)*rhov_T + EoS.d2ud2T(fv3)*rhov_T^2 + EoS.dudT(fv3)*rhov_TT;
    sl_TT := EoS.d2sT2d(fl3) + 2*EoS.d2sTd(fl3)*rhol_T + EoS.d2sd2T(fl3)*rhol_T^2 + EoS.dsdT(fl3)*rhol_TT;
    sv_TT := EoS.d2sT2d(fv3) + 2*EoS.d2sTd(fv3)*rhov_T + EoS.d2sd2T(fv3)*rhov_T^2 + EoS.dsdT(fv3)*rhov_TT;
    vl_TT := -rhol_TT*vl^2 + 2*rhol_T^2*vl^3;
    vv_TT := -rhov_TT*vv^2 + 2*rhov_T^2*vv^3;
    x := (1/d - vl)/(vv - vl);
    x_d := -1/(d^2*(vv - vl));
    x_T := (-vl_T*(vv - vl) - (1/d - vl)*(vv_T - vl_T))/(vv - vl)^2;
    // two phase: a = sqrt(1/ddps), ddps = (dd/dp)@s = -d^2*G/P1, with
    // G = dv'/dT + (x*ds"/dT + (1-x)*ds'/dT)*(v"-v')/(s'-s") + x*(dv"/dT - dv'/dT) (as in velocityOfSound_cache)
    G := vl_T + (x*sv_T + (1 - x)*sl_T)*(vv - vl)/(sl - sv) + x*(vv_T - vl_T);
    G_x := (sv_T - sl_T)*(vv - vl)/(sl - sv) + (vv_T - vl_T);
    G_T := vl_TT + (x*sv_TT + (1 - x)*sl_TT)*(vv - vl)/(sl - sv)
         + (x*sv_T + (1 - x)*sl_T)*((vv_T - vl_T)*(sl - sv) - (vv - vl)*(sl_T - sv_T))/(sl - sv)^2
         + x*(vv_TT - vl_TT);
    ddps := -d^2*G/P1;
    ddps_t := (-2*d*G/P1 - d^2*G_x*x_d/P1)*d_der + (-d^2*(G_T + G_x*x_T)/P1 + d^2*G*P2/P1^2)*T_der;
    a_der := -0.5*ddps_t/ddps^1.5;
  else
    assert(false, "velocityOfSound_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.velocityOfSound_cache\">velocityOfSound_cache</a>,
used for analytic Jacobians. In the two-phase region the velocity of sound follows from (dd/dp)@s, which depends on d and on the saturation properties and their first derivatives along the saturation line; its derivative needs the second derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end velocityOfSound_cache_der;
