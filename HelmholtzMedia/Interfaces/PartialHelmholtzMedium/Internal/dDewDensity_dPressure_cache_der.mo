within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function dDewDensity_dPressure_cache_der "time derivative of dDewDensity_dPressure_cache"
  extends Modelica.Icons.Function;
  input Temperature T "saturation temperature, sat.Tsat";
  input SaturationProperties sat;
  input EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at the saturated state (not used)";
  input Real T_der "time derivative of T";
  output Real ddvdp_der "time derivative of dDewDensity_dPressure";
protected
  EoS.HelmholtzDerivs fl3=EoS.setHelmholtzDerivs(d=sat.liq.d, T=T, order=3) "Helmholtz derivatives up to third order at the saturated liquid";
  EoS.HelmholtzDerivs fv3=EoS.setHelmholtzDerivs(d=sat.vap.d, T=T, order=3) "Helmholtz derivatives up to third order at the saturated vapour";
  Real vl=1/sat.liq.d "specific volume of saturated liquid";
  Real vv=1/sat.vap.d "specific volume of saturated vapour";
  Real sl=sat.liq.s "specific entropy of saturated liquid";
  Real sv=sat.vap.s "specific entropy of saturated vapour";
  // Clausius-Clapeyron
  Real P1=(sv - sl)/(vv - vl) "dp_sat/dT along the saturation line";
  Real rhol_T=(P1 - EoS.dpTd(fl3))/EoS.dpdT(fl3) "d(rho_liq)/dT along saturation";
  Real rhov_T=(P1 - EoS.dpTd(fv3))/EoS.dpdT(fv3) "d(rho_vap)/dT along saturation";
  Real sl_T=EoS.dsTd(fl3) + EoS.dsdT(fl3)*rhol_T "d(s_liq)/dT along saturation";
  Real sv_T=EoS.dsTd(fv3) + EoS.dsdT(fv3)*rhov_T "d(s_vap)/dT along saturation";
  Real vl_T=-rhol_T*vl^2 "d(v_liq)/dT along saturation";
  Real vv_T=-rhov_T*vv^2 "d(v_vap)/dT along saturation";
  Real P2=((sv_T - sl_T)*(vv - vl) - (sv - sl)*(vv_T - vl_T))/(vv - vl)^2 "d2p_sat/dT2 along the saturation line";
  Real rhov_TT=(P2 - EoS.d2pT2d(fv3) - EoS.d2pTd(fv3)*rhov_T)/EoS.dpdT(fv3)
      - (P1 - EoS.dpTd(fv3))*(EoS.d2pTd(fv3) + EoS.d2pd2T(fv3)*rhov_T)/EoS.dpdT(fv3)^2;
algorithm
  // the value is d(rhov)/dp along saturation = rhov_T/P1
  ddvdp_der := (rhov_TT*P1 - rhov_T*P2)/P1^2*T_der;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.dDewDensity_dPressure_cache\">dDewDensity_dPressure_cache</a>,
used for analytic Jacobians. Along the saturation line the value depends on T only; its derivative needs the second derivatives along the saturation line.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end dDewDensity_dPressure_cache_der;
