within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derp_h_cache_der "time derivative of density_derp_h_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real ddph_der "time derivative of (dd/dp)@h";
protected
  EoS.HelmholtzDerivs f3 "Helmholtz derivatives up to third order at (d, T)";
  Real A "(dp/dd)@T";
  Real B "(dp/dT)@d";
  Real C "(dh/dd)@T";
  Real D "(dh/dT)@d";
  Real A_t "time derivative of A";
  Real B_t "time derivative of B";
  Real C_t "time derivative of C";
  Real D_t "time derivative of D";
  Real F;
  Real F_t;
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
  EoS.HelmholtzDerivs fl3 "Helmholtz derivatives up to third order at the saturated liquid";
  EoS.HelmholtzDerivs fv3 "Helmholtz derivatives up to third order at the saturated vapour";
  Real P2 "d2p_sat/dT2 along the saturation line";
  Real sl_T "d(s_liq)/dT along saturation";
  Real sv_T "d(s_vap)/dT along saturation";
  Real rhol_TT "d2(rho_liq)/dT2 along saturation";
  Real rhov_TT "d2(rho_vap)/dT2 along saturation";
  Real hl_TT "d2(h_liq)/dT2 along saturation";
  Real hv_TT "d2(h_vap)/dT2 along saturation";
  Real vl_TT "d2(v_liq)/dT2 along saturation";
  Real vv_TT "d2(v_vap)/dT2 along saturation";
  Real x "vapour quality";
  Real x_d "(dx/dd)@T";
  Real x_T "(dx/dT)@d";
  Real G "(dv/dp)@h*dp_sat/dT";
  Real G_x "(dG/dx)@T";
  Real G_T "(dG/dT)@x";
  Real dd "partial derivative of ddph by d at constant T";
  Real dT "partial derivative of ddph by T at constant d";
algorithm
  if phase==1 then
    f3 := EoS.setHelmholtzDerivs(d=d, T=T, order=3);
    A := EoS.dpdT(f3);
    B := EoS.dpTd(f3);
    C := EoS.dhdT(f3);
    D := EoS.dhTd(f3);
    A_t := EoS.d2pd2T(f3)*d_der + EoS.d2pTd(f3)*T_der;
    B_t := EoS.d2pTd(f3)*d_der + EoS.d2pT2d(f3)*T_der;
    C_t := EoS.d2hd2T(f3)*d_der + EoS.d2hTd(f3)*T_der;
    D_t := EoS.d2hTd(f3)*d_der + EoS.d2hT2d(f3)*T_der;
    // ddph = 1/F with F = A - B*C/D
    F := A - B*C/D;
    F_t := A_t - (B_t*C + B*C_t)/D + B*C*D_t/D^2;
    ddph_der := -F_t/F^2;
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
    // second derivatives along the saturation line
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
    // ddph = -d^2*G/P1, with G = dv'/dT + (x*dh"/dT + (1-x)*dh'/dT)*(v"-v')/(h'-h") + x*(dv"/dT - dv'/dT)
    // (the same expression as in density_derp_h_cache, written with derivatives along saturation by T)
    x := (1/d - vl)/(vv - vl);
    x_d := -1/(d^2*(vv - vl));
    x_T := (-vl_T*(vv - vl) - (1/d - vl)*(vv_T - vl_T))/(vv - vl)^2;
    G := vl_T + (x*hv_T + (1 - x)*hl_T)*(vv - vl)/(hl - hv) + x*(vv_T - vl_T);
    G_x := (hv_T - hl_T)*(vv - vl)/(hl - hv) + (vv_T - vl_T);
    G_T := vl_TT + (x*hv_TT + (1 - x)*hl_TT)*(vv - vl)/(hl - hv)
         + (x*hv_T + (1 - x)*hl_T)*((vv_T - vl_T)*(hl - hv) - (vv - vl)*(hl_T - hv_T))/(hl - hv)^2
         + x*(vv_TT - vl_TT);
    dd := -2*d*G/P1 - d^2*G_x*x_d/P1;
    dT := -d^2*(G_T + G_x*x_T)/P1 + d^2*G*P2/P1^2;
    ddph_der := dd*d_der + dT*T_der;
  else
    assert(false, "density_derp_h_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.density_derp_h_cache\">density_derp_h_cache</a>,
used for analytic Jacobians. In the two-phase region (dd/dp)@h depends on d and, through the saturation properties and their first derivatives along the saturation line, on T; its time derivative needs the second derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end density_derp_h_cache_der;
