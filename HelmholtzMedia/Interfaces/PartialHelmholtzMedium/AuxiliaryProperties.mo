within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
record AuxiliaryProperties
  "intermediate results shared by the property functions of one state, see setAuxiliaryProperties"
  extends Modelica.Icons.Record;
  EoS.HelmholtzDerivs f "derivatives of the Helmholtz energy at the state (single phase)";
  SaturationProperties sat "saturation properties at the temperature of the state (two phase)";
  EoS.HelmholtzDerivs fl "derivatives of the Helmholtz energy at the saturated liquid density (two phase)";
  EoS.HelmholtzDerivs fv "derivatives of the Helmholtz energy at the saturated vapour density (two phase)";
end AuxiliaryProperties;
