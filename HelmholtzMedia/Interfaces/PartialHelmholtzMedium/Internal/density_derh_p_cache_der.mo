within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function density_derh_p_cache_der "time derivative of density_derh_p_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real ddhp_der "time derivative of (dd/dh)@p";
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
  Real E;
  Real E_t;
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
  Real Q;
  Real Q_T;
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
    // ddhp = 1/E with E = C - D*A/B
    E := C - D*A/B;
    E_t := C_t - (D_t*A + D*A_t)/B + D*A*B_t/B^2;
    ddhp_der := -E_t/E^2;
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
    // ddhp = -d^2*Q with Q = (v'-v")/(h'-h"), saturation properties depend on T only
    Q := (vl - vv)/(hl - hv);
    Q_T := ((vl_T - vv_T)*(hl - hv) - (vl - vv)*(hl_T - hv_T))/(hl - hv)^2;
    ddhp_der := -2*d*d_der*Q - d^2*Q_T*T_der;
  else
    assert(false, "density_derh_p_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.density_derh_p_cache\">density_derh_p_cache</a>,
used for analytic Jacobians. In the two-phase region (dd/dh)@p = -d^2*(v_liq - v_vap)/(h_liq - h_vap) depends on d and on the saturation properties at T.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end density_derh_p_cache_der;
