within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_dilute
  "Return thermal conductivity dilute contribution"
  // depends on dynamicViscosity, specificHeatCapacityCp, specificHeatCapacityCv and dpdd=1/dddp
  input ThermodynamicState state;
  output ThermalConductivity lambda_0;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_residual=fluidConstants[1].molarMass/
      thermalConductivityCoefficients.reducingMolarVolume_residual;
  Real delta "reduced density";

  Temperature T_red_0=thermalConductivityCoefficients.reducingTemperature_0;
  Real tau "reduced temperature";

  // coeffs for dilute contribution
  Integer nDilute_num = size(thermalConductivityCoefficients.lambda_0_num_coeffs,1);
  Real[nDilute_num, 2] A_num= thermalConductivityCoefficients.lambda_0_num_coeffs;
  Integer nDilute_den = size(thermalConductivityCoefficients.lambda_0_den_coeffs,1);
  Real[nDilute_den, 2] A_den= thermalConductivityCoefficients.lambda_0_den_coeffs;
  Real denom= 1;
  Real lambda_red_0=thermalConductivityCoefficients.reducingThermalConductivity_0;

  // interim variables for flagged algorithms
  constant Real eps = Modelica.Constants.eps;
  EoS.HelmholtzDerivs f;
  Real wm = 1000*MM; // RefProp uses g per mol
  Real cp0=1;  // ideal gas cp
  Real eta_0=1; // dilute contribution only

algorithm
  Modelica.Utilities.Streams.print("thermalConductivity_dilute: d = " + String(state.d) + " and T = " + String(state.T));
  tau := state.T/T_red_0;

  // numerator terms, check exponent
  if (A_num[nDilute_num,2]>-90) then
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num);
  else
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num-1);
    f:=EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);
    cp0:=f.R*(1 - f.tau*f.tau*f.itt);
    eta_0:=dynamicViscosity_dilute(state);
    Modelica.Utilities.Streams.print("flagged algorithms, cp0=" + String(cp0) + " and eta_0=" + String(eta_0));
    if (abs(A_num[nDilute_num,2]+99)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 99: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+98)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 98: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+97)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 97: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+96)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 96");
      cp0 := cp0/f.R-2.5;
      lambda_0 := (lambda_0*cp0+15.0/4.0)*f.R*eta_0/wm;
    end if;
  end if;

  if (nDilute_den>1) then
    // Modelica.Utilities.Streams.print("denominator term");
    denom := sum(A_den[i, 1]*tau^A_den[i, 2] for i in 1:nDilute_den);
  else
    denom := 1;
  end if;

  lambda_0 := lambda_0/denom*lambda_red_0;
  Modelica.Utilities.Streams.print(" lambda_0 = " + String(lambda_0));

end thermalConductivity_dilute;
