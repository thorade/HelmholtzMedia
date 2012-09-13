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
  constant Real eps = Modelica.Constants.eps;
  Real lambda_red_0=thermalConductivityCoefficients.reducingThermalConductivity_0;

  EoS.HelmholtzDerivs f;
  SpecificHeatCapacity cp0=1;
  DynamicViscosity eta_0=1;

algorithm
  // dilute gas contribution
  tau := state.T/T_red_0;

  // numerator terms, check exponent
  if (A_num[nDilute_num,2]>-90) then
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num);
  else
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num-1);
    f:=EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);
    cp0:=f.R*(1 - f.tau*f.tau*f.itt);
    eta_0:=dynamicViscosity_dilute(state);
    if (abs(A_num[nDilute_num,2]+99)<eps) then
      Modelica.Utilities.Streams.print("flag 99: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+98)<eps) then
      Modelica.Utilities.Streams.print("flag 98: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+97)<eps) then
      Modelica.Utilities.Streams.print("flag 97: not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+96)<eps) then
      Modelica.Utilities.Streams.print("flag 96");
      //lambda_0 := (lambda_0*(cp0/f.R-2.5)+15.0/4.0)*f.R*eta_0/MM;
      lambda_0 := lambda_0*eta_0/MM*(cp0-2.5*f.R) + 15.0/4.0*f.R*eta_0/MM;
    end if;
  end if;

  if (nDilute_den>1) then
    // Modelica.Utilities.Streams.print("denominator term");
    denom := sum(A_den[i, 1]*tau^A_den[i, 2] for i in 1:nDilute_den);
  else
    denom := 1;
  end if;

  lambda_0 := lambda_0/denom*lambda_red_0;

   // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print("      cp0 = " + String(cp0));
  Modelica.Utilities.Streams.print("    eta_0 = " + String(eta_0));
  Modelica.Utilities.Streams.print(" lambda_0 = " + String(lambda_0));
  Modelica.Utilities.Streams.print("===========================================");

end thermalConductivity_dilute;
