within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function specificHeatCapacityCp_cache_der "time derivative of specificHeatCapacityCp_cache"
  extends Modelica.Icons.Function;
  input Density d "density of the state";
  input Temperature T "temperature of the state";
  input FixedPhase phase "phase of the state: 1 single phase, 2 two-phase";
  input HelmholtzCache cache "see setHelmholtzCache(d, T, phase)";
  input Real d_der "time derivative of d";
  input Real T_der "time derivative of T";
  output Real cp_der "time derivative of the isobaric specific heat capacity";
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
    // cp = D - C*B/A
    cp_der := D_t - (C_t*B + C*B_t)/A + C*B*A_t/A^2;
  elseif phase==2 then
    cp_der := 0 "cp is infinite in the two-phase region";
  else
    assert(false, "specificHeatCapacityCp_cache_der error: phase is neither 1 nor 2");
  end if;
  annotation (Documentation(info="<html>
<p>
Time derivative of <a href=\"modelica://HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal.specificHeatCapacityCp_cache\">specificHeatCapacityCp_cache</a>,
used for analytic Jacobians. In the two-phase region cp is infinite, its time derivative is set to zero.
The Helmholtz derivatives up to third order are computed here, the cache (second order) is not differentiated.
</p>
</html>"));
end specificHeatCapacityCp_cache_der;
