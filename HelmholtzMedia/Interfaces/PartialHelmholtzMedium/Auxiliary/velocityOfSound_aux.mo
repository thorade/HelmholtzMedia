within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function velocityOfSound_aux "returns the speed or velocity of sound (computed from the shared auxiliary properties, see velocityOfSound)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output VelocityOfSound a;
protected
  MassFraction x "vapour quality";
  DerTemperatureByPressure dTp;
  DerEntropyByPressure dsp_liq;
  DerEntropyByPressure dsp_vap;
  DerDensityByPressure ddp_liq;
  DerDensityByPressure ddp_vap;
  DerVolumeByPressure dvp_liq;
  DerVolumeByPressure dvp_vap;
  DerVolumeByPressure dvps;
  DerFractionByPressure dxps;
algorithm
  if state.phase==1 then
    a := sqrt(EoS.dpdT(aux.f)-EoS.dpTd(aux.f)*EoS.dsdT(aux.f)/EoS.dsTd(aux.f));
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=aux.sat);
    dTp := saturationTemperature_derp_sat(sat=aux.sat);
    dsp_liq := EoS.dsdT(aux.fl)/EoS.dpdT(aux.fl) + (EoS.dsTd(aux.fl)-EoS.dsdT(aux.fl)*EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl))*dTp;
    dsp_vap := EoS.dsdT(aux.fv)/EoS.dpdT(aux.fv) + (EoS.dsTd(aux.fv)-EoS.dsdT(aux.fv)*EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv))*dTp;
    ddp_liq := 1.0/EoS.dpdT(aux.fl) - EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl)*dTp;
    ddp_vap := 1.0/EoS.dpdT(aux.fv) - EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv)*dTp;
    dvp_liq := -1.0/aux.sat.liq.d^2 * ddp_liq;
    dvp_vap := -1.0/aux.sat.vap.d^2 * ddp_vap;
    dxps :=(x*dsp_vap + (1 - x)*dsp_liq)/(aux.sat.liq.s - aux.sat.vap.s);
    dvps := dvp_liq + dxps*(1.0/aux.sat.vap.d-1.0/aux.sat.liq.d) + x*(dvp_vap-dvp_liq);
    a := sqrt(-1.0/(state.d*state.d*dvps));
  end if;
end velocityOfSound_aux;
