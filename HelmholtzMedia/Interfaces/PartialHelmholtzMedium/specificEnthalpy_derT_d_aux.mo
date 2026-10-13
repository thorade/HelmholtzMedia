within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d_aux "returns enthalpy derivative (dh/dT)@d=const (computed from the shared auxiliary properties, see specificEnthalpy_derT_d)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output DerEnthalpyByTemperature dhTd;
protected
  MassFraction x "vapour quality";
  DerPressureByTemperature dpT;
  DerEnthalpyByTemperature dhT_liq;
  DerEnthalpyByTemperature dhT_vap;
  DerDensityByTemperature ddT_liq;
  DerDensityByTemperature ddT_vap;
  DerVolumeByTemperature dvT_liq;
  DerVolumeByTemperature dvT_vap;
  DerFractionByTemperature dxTv;
algorithm
  if state.phase==1 then
    dhTd := EoS.dhTd(aux.f);
  elseif state.phase==2 then
    x := (1/state.d - 1/aux.sat.liq.d)/(1/aux.sat.vap.d - 1/aux.sat.liq.d);
    dpT := (aux.sat.vap.s-aux.sat.liq.s)/(1.0/aux.sat.vap.d-1.0/aux.sat.liq.d);
    dhT_liq := EoS.dhTd(aux.fl)-EoS.dhdT(aux.fl)*EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl) + EoS.dhdT(aux.fl)/EoS.dpdT(aux.fl)*dpT;
    dhT_vap := EoS.dhTd(aux.fv)-EoS.dhdT(aux.fv)*EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv) + EoS.dhdT(aux.fv)/EoS.dpdT(aux.fv)*dpT;
    ddT_liq := -EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl) + 1.0/EoS.dpdT(aux.fl)*dpT;
    ddT_vap := -EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv) + 1.0/EoS.dpdT(aux.fv)*dpT;
    dvT_liq := -1/aux.sat.liq.d^2 * ddT_liq;
    dvT_vap := -1/aux.sat.vap.d^2 * ddT_vap;
    dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/aux.sat.liq.d - 1/aux.sat.vap.d);
    dhTd := dhT_liq + dxTv*(aux.sat.vap.h-aux.sat.liq.h) + x*(dhT_vap-dhT_liq);
  else
    assert(false, "specificEnthalpy_derT_d error: phase is neither 1 nor 2");
  end if;
end specificEnthalpy_derT_d_aux;
