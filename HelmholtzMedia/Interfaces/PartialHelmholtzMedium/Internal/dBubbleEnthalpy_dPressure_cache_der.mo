within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dBubbleEnthalpy_dPressure_cache_der "time derivative of dBubbleEnthalpy_dPressure_cache"
  extends Modelica.Icons.Function;
  input Temperature T "saturation temperature, sat.Tsat";
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at the saturated state (not used)";
  input Real T_der "time derivative of T";
  output Real dhldp_der "time derivative of dBubbleEnthalpy_dPressure";
protected
  EoS.HelmholtzDerivs fl3=EoS.setHelmholtzDerivs(d=sat.liq.d, T=T, order=3) "Helmholtz derivatives up to third order at the saturated liquid";
  EoS.HelmholtzDerivs fv3=EoS.setHelmholtzDerivs(d=sat.vap.d, T=T, order=3) "Helmholtz derivatives up to third order at the saturated vapour";
  Real vl=1/sat.liq.d "specific volume of saturated liquid";
  Real vv=1/sat.vap.d "specific volume of saturated vapour";
  Real hl=sat.liq.h "specific enthalpy of saturated liquid";
  Real hv=sat.vap.h "specific enthalpy of saturated vapour";
  Real ul=sat.liq.u "specific inner energy of saturated liquid";
  Real uv=sat.vap.u "specific inner energy of saturated vapour";
  Real sl=sat.liq.s "specific entropy of saturated liquid";
  Real sv=sat.vap.s "specific entropy of saturated vapour";
  // Clausius-Clapeyron
  Real P1=(sv - sl)/(vv - vl) "dp_sat/dT along the saturation line";
  Real rhol_T=(P1 - EoS.dpTd(fl3))/EoS.dpdT(fl3) "d(rho_liq)/dT along saturation";
  Real rhov_T=(P1 - EoS.dpTd(fv3))/EoS.dpdT(fv3) "d(rho_vap)/dT along saturation";
  Real hl_T=EoS.dhTd(fl3) + EoS.dhdT(fl3)*rhol_T "d(h_liq)/dT along saturation";
  Real hv_T=EoS.dhTd(fv3) + EoS.dhdT(fv3)*rhov_T "d(h_vap)/dT along saturation";
  Real ul_T=EoS.duTd(fl3) + EoS.dudT(fl3)*rhol_T "d(u_liq)/dT along saturation";
  Real uv_T=EoS.duTd(fv3) + EoS.dudT(fv3)*rhov_T "d(u_vap)/dT along saturation";
  Real sl_T=EoS.dsTd(fl3) + EoS.dsdT(fl3)*rhol_T "d(s_liq)/dT along saturation";
  Real sv_T=EoS.dsTd(fv3) + EoS.dsdT(fv3)*rhov_T "d(s_vap)/dT along saturation";
  Real vl_T=-rhol_T*vl^2 "d(v_liq)/dT along saturation";
  Real vv_T=-rhov_T*vv^2 "d(v_vap)/dT along saturation";
  Real P2=((sv_T - sl_T)*(vv - vl) - (sv - sl)*(vv_T - vl_T))/(vv - vl)^2 "d2p_sat/dT2 along the saturation line";
  Real rhol_TT=(P2 - EoS.d2pT2d(fl3) - EoS.d2pTd(fl3)*rhol_T)/EoS.dpdT(fl3)
      - (P1 - EoS.dpTd(fl3))*(EoS.d2pTd(fl3) + EoS.d2pd2T(fl3)*rhol_T)/EoS.dpdT(fl3)^2;
  Real rhov_TT=(P2 - EoS.d2pT2d(fv3) - EoS.d2pTd(fv3)*rhov_T)/EoS.dpdT(fv3)
      - (P1 - EoS.dpTd(fv3))*(EoS.d2pTd(fv3) + EoS.d2pd2T(fv3)*rhov_T)/EoS.dpdT(fv3)^2;
  Real hl_TT=EoS.d2hT2d(fl3) + 2*EoS.d2hTd(fl3)*rhol_T + EoS.d2hd2T(fl3)*rhol_T^2 + EoS.dhdT(fl3)*rhol_TT;
  Real hv_TT=EoS.d2hT2d(fv3) + 2*EoS.d2hTd(fv3)*rhov_T + EoS.d2hd2T(fv3)*rhov_T^2 + EoS.dhdT(fv3)*rhov_TT;
  Real ul_TT=EoS.d2uT2d(fl3) + 2*EoS.d2uTd(fl3)*rhol_T + EoS.d2ud2T(fl3)*rhol_T^2 + EoS.dudT(fl3)*rhol_TT;
  Real uv_TT=EoS.d2uT2d(fv3) + 2*EoS.d2uTd(fv3)*rhov_T + EoS.d2ud2T(fv3)*rhov_T^2 + EoS.dudT(fv3)*rhov_TT;
  Real sl_TT=EoS.d2sT2d(fl3) + 2*EoS.d2sTd(fl3)*rhol_T + EoS.d2sd2T(fl3)*rhol_T^2 + EoS.dsdT(fl3)*rhol_TT;
  Real sv_TT=EoS.d2sT2d(fv3) + 2*EoS.d2sTd(fv3)*rhov_T + EoS.d2sd2T(fv3)*rhov_T^2 + EoS.dsdT(fv3)*rhov_TT;
  Real vl_TT=-rhol_TT*vl^2 + 2*rhol_T^2*vl^3;
  Real vv_TT=-rhov_TT*vv^2 + 2*rhov_T^2*vv^3;
algorithm
  // the value is d(hl)/dp along saturation = hl_T/P1
  dhldp_der := (hl_TT*P1 - hl_T*P2)/P1^2*T_der;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.dBubbleEnthalpy_dPressure_cache\">dBubbleEnthalpy_dPressure_cache</a>,
used for analytic Jacobians. Along the saturation line the value depends on T only; its derivative needs the second derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end dBubbleEnthalpy_dPressure_cache_der;
