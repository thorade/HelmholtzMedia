within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function setSat_T_der "time derivative of setSat_T (derivative of the saturation record)"
  extends Modelica.Icons.Function;
  input Temperature T;
  input Real T_der "time derivative of T";
  output SaturationProperties sat_der "time derivative of the saturation state";
protected
  SaturationProperties sat "the saturation state (computed again)";
  EoS.HelmholtzDerivs fl "Helmholtz derivatives at the saturated liquid";
  EoS.HelmholtzDerivs fv "Helmholtz derivatives at the saturated vapour";
  Real Tsat_der "time derivative of the saturation temperature";
  Real dl_der "time derivative of the liquid density";
  Real dv_der "time derivative of the vapour density";
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
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
algorithm
  if (T_der == 0) then
    // zero seed (e.g. in the analytic Jacobian of a nonlinear system): skip the evaluation of the saturation state
    sat_der.psat := 0;
    sat_der.Tsat := 0;
    sat_der.liq := ThermodynamicState(phase=0, d=0, T=0, p=0, h=0, u=0, s=0);
    sat_der.vap := ThermodynamicState(phase=0, d=0, T=0, p=0, h=0, u=0, s=0);
    return;
  end if;
  sat := setSat_T(T=T);
  fl := EoS.setHelmholtzDerivs(d=sat.liq.d, T=sat.Tsat, order=2);
  fv := EoS.setHelmholtzDerivs(d=sat.vap.d, T=sat.Tsat, order=2);
  if ((T>=T_trip) and (T<T_crit*(1-1e-7))) then
    vl := 1/sat.liq.d;
    vv := 1/sat.vap.d;
    hl := sat.liq.h;
    hv := sat.vap.h;
    ul := sat.liq.u;
    uv := sat.vap.u;
    sl := sat.liq.s;
    sv := sat.vap.s;
    P1 := (sv - sl)/(vv - vl) "Clausius-Clapeyron";
    rhol_T := (P1 - EoS.dpTd(fl))/EoS.dpdT(fl);
    rhov_T := (P1 - EoS.dpTd(fv))/EoS.dpdT(fv);
    hl_T := EoS.dhTd(fl) + EoS.dhdT(fl)*rhol_T;
    hv_T := EoS.dhTd(fv) + EoS.dhdT(fv)*rhov_T;
    ul_T := EoS.duTd(fl) + EoS.dudT(fl)*rhol_T;
    uv_T := EoS.duTd(fv) + EoS.dudT(fv)*rhov_T;
    sl_T := EoS.dsTd(fl) + EoS.dsdT(fl)*rhol_T;
    sv_T := EoS.dsTd(fv) + EoS.dsdT(fv)*rhov_T;
    vl_T := -rhol_T*vl^2;
    vv_T := -rhov_T*vv^2;
    P2 := ((sv_T - sl_T)*(vv - vl) - (sv - sl)*(vv_T - vl_T))/(vv - vl)^2;
    dl_der := rhol_T*T_der;
    dv_der := rhov_T*T_der;
  elseif (T>=T_crit*(1-1e-7)) then
    // critical isochore
    dl_der := 0;
    dv_der := 0;
  else
    // below the triple point: d_liq ~ 1/T, d_vap ~ T (see setSat_T)
    dl_der := -sat.liq.d/sat.Tsat*T_der;
    dv_der := sat.vap.d/sat.Tsat*T_der;
  end if;
  sat_der.Tsat := T_der;
  sat_der.liq := stateDerivative_dT(f=fl, d_der=dl_der, T_der=T_der);
  sat_der.vap := stateDerivative_dT(f=fv, d_der=dv_der, T_der=T_der);
  sat_der.psat := (sat_der.liq.p + sat_der.vap.p)/2 "as in setSat_T";
  annotation (Documentation(info="<html>
<p>
Derivative annotation of setSat_T. The output has the same record type as the saturation properties, with the time
derivatives of the Real components (phase=0), see setState_phX_der. Along the saturation line the densities change
with d(rho)/dT = (dp_sat/dT - (dp/dT)@d)/(dp/dd)@T, with dp_sat/dT from Clausius-Clapeyron.
</p>
</html>"));
end setSat_T_der;
