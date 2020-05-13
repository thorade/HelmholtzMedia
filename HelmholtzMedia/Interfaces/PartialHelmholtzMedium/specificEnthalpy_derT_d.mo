within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function specificEnthalpy_derT_d "returns enthalpy derivative (dh/dT)@d=const"
  input ThermodynamicState state;
  output DerEnthalpyByTemperature dhTd;

protected
  EoS.HelmholtzDerivs f;

  SaturationProperties sat;
  MassFraction x "vapour quality";
  DerPressureByTemperature dpT;
  EoS.HelmholtzDerivs fl;
  EoS.HelmholtzDerivs fv;

  DerEnthalpyByTemperature dhT_liq;
  DerEnthalpyByTemperature dhT_vap;
  DerDensityByTemperature ddT_liq;
  DerDensityByTemperature ddT_vap;
  DerVolumeByTemperature dvT_liq;
  DerVolumeByTemperature dvT_vap;
  DerFractionByTemperature dxTv;

algorithm
  if state.phase==1 then
    f := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=state.phase);
    dhTd := EoS.dhTd(f);

  elseif state.phase==2 then
    sat:=setSat_T(T=state.T);
    x := (1/state.d - 1/sat.liq.d)/(1/sat.vap.d - 1/sat.liq.d);
    dpT := (sat.vap.s-sat.liq.s)/(1.0/sat.vap.d-1.0/sat.liq.d);

    fl := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.liq.d, phase=1);
    fv := EoS.setHelmholtzDerivsSecond(T=state.T, d=sat.vap.d, phase=1);

    dhT_liq := EoS.dhTd(fl)-EoS.dhdT(fl)*EoS.dpTd(fl)/EoS.dpdT(fl) + EoS.dhdT(fl)/EoS.dpdT(fl)*dpT;
    dhT_vap := EoS.dhTd(fv)-EoS.dhdT(fv)*EoS.dpTd(fv)/EoS.dpdT(fv) + EoS.dhdT(fv)/EoS.dpdT(fv)*dpT;
    ddT_liq := -EoS.dpTd(fl)/EoS.dpdT(fl) + 1.0/EoS.dpdT(fl)*dpT;
    ddT_vap := -EoS.dpTd(fv)/EoS.dpdT(fv) + 1.0/EoS.dpdT(fv)*dpT;
    dvT_liq := -1/sat.liq.d^2 * ddT_liq;
    dvT_vap := -1/sat.vap.d^2 * ddT_vap;
    dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/sat.liq.d - 1/sat.vap.d);

    dhTd := dhT_liq + dxTv*(sat.vap.h-sat.liq.h) + x*(dhT_vap-dhT_liq);
  end if;

end specificEnthalpy_derT_d;
