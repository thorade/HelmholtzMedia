within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function isothermalCompressibility_cache_der "time derivative of isothermalCompressibility_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real kappa_der "time derivative of the isothermal compressibility";
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
    kappa_der := -(d_der*A + d*A_t)/(d*A)^2 "kappa = 1/(d*A)";
  elseif phase==2 then
    kappa_der := 0 "infinite in the two-phase region";
  else
    assert(false, "isothermalCompressibility_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.isothermalCompressibility_cache\">isothermalCompressibility_cache</a>,
used for analytic Jacobians. In the two-phase region the value is infinite, its derivative is set to zero.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end isothermalCompressibility_cache_der;
