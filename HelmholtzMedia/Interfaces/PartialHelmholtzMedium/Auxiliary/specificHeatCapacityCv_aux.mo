within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Auxiliary;
function specificHeatCapacityCv_aux "returns the isochoric specific heat capcacity (computed from the shared auxiliary properties, see specificHeatCapacityCv)"
  input ThermodynamicState state;
  input AuxiliaryProperties aux;
  output SpecificHeatCapacity cv;
protected
  MassFraction x "vapour quality";
  DerPressureByTemperature dpT;
  DerEnergyByTemperature duT_liq;
  DerEnergyByTemperature duT_vap;
  DerDensityByTemperature ddT_liq;
  DerDensityByTemperature ddT_vap;
  DerVolumeByTemperature dvT_liq;
  DerVolumeByTemperature dvT_vap;
  DerFractionByTemperature dxTv;
algorithm
  if state.phase==1 then
    cv := EoS.duTd(aux.f);
  elseif state.phase==2 then
    x := vapourQuality_sat(state=state, sat=aux.sat);
    dpT := saturationPressure_derT_sat(sat=aux.sat);
    duT_liq := EoS.duTd(aux.fl)-EoS.dudT(aux.fl)*EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl) + EoS.dudT(aux.fl)/EoS.dpdT(aux.fl)*dpT;
    duT_vap := EoS.duTd(aux.fv)-EoS.dudT(aux.fv)*EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv) + EoS.dudT(aux.fv)/EoS.dpdT(aux.fv)*dpT;
    ddT_liq := -EoS.dpTd(aux.fl)/EoS.dpdT(aux.fl) + 1.0/EoS.dpdT(aux.fl)*dpT;
    ddT_vap := -EoS.dpTd(aux.fv)/EoS.dpdT(aux.fv) + 1.0/EoS.dpdT(aux.fv)*dpT;
    dvT_liq := -1/aux.sat.liq.d^2 * ddT_liq;
    dvT_vap := -1/aux.sat.vap.d^2 * ddT_vap;
    dxTv :=(x*dvT_vap + (1 - x)*dvT_liq)/(1/aux.sat.liq.d - 1/aux.sat.vap.d);
    cv := duT_liq + dxTv*(aux.sat.vap.u-aux.sat.liq.u) + x*(duT_vap-duT_liq);
  end if;
end specificHeatCapacityCv_aux;
