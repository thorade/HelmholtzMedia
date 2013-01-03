within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Ancillary;
function meltingPressure_T
  "ancillary function: calculate melting pressure for a given Temperature"
  input Modelica.SIunits.Temperature T;
  output Modelica.SIunits.AbsolutePressure p_melt;

protected
  PressureMeltingModel pressureMeltingModel=ancillaryCoefficients.pressureMeltingModel;
  Temperature T_reducing=ancillaryCoefficients.T_reducing;
  AbsolutePressure p_reducing=ancillaryCoefficients.p_reducing;

  Integer nPressureMelting1=size(ancillaryCoefficients.pressureMelting1, 1);
  Real[nPressureMelting1,2] n1=ancillaryCoefficients.pressureMelting1;

  Integer nPressureMelting2=size(ancillaryCoefficients.pressureMelting2, 1);
  Real[nPressureMelting2,2] n2=ancillaryCoefficients.pressureMelting2;

  Integer nPressureMelting3=size(ancillaryCoefficients.pressureMelting3, 1);
  Real[nPressureMelting3,2] n3=ancillaryCoefficients.pressureMelting3;

  constant Real Tr=max((T/T_reducing), Modelica.Constants.small)
    "reduced temperature";
  Real pr "reduced pressure";

algorithm
  pr := 0;
  pr := sum(n1[i,1]*Tr^n1[i,2] for i in 1:nPressureMelting1)
      + sum(n2[i,1]*(Tr-1)^n2[i,2] for i in 1:nPressureMelting2)
      + sum(n3[i,1]*log(Tr)^n3[i,2] for i in 1:nPressureMelting3);

  if (pressureMeltingModel == PressureMeltingModel.ML1) then
    p_melt := p_reducing*pr;
  elseif (pressureMeltingModel == PressureMeltingModel.ML2) then
    p_melt := p_reducing*exp(pr);
  end if;

  if (pr<1e-6) then
    p_melt := fluidLimits.PMAX;
  end if;

end meltingPressure_T;
