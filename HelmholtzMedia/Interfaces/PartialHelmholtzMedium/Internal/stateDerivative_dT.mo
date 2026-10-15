within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
function stateDerivative_dT "time derivative of a single-phase state, given the time derivatives of d and T"
  extends Modelica.Icons.Function;
  input EoS.HelmholtzDerivs f "Helmholtz derivatives (at least second order) at the state";
  input Real d_der "time derivative of density";
  input Real T_der "time derivative of temperature";
  output ThermodynamicState state_der "time derivative of the state (phase=0)";
algorithm
  state_der := ThermodynamicState(
    phase=0,
    T=T_der,
    p=EoS.dpdT(f)*d_der + EoS.dpTd(f)*T_der,
    d=d_der,
    u=EoS.dudT(f)*d_der + EoS.duTd(f)*T_der,
    h=EoS.dhdT(f)*d_der + EoS.dhTd(f)*T_der,
    s=EoS.dsdT(f)*d_der + EoS.dsTd(f)*T_der);
  annotation (Inline=true);
end stateDerivative_dT;
