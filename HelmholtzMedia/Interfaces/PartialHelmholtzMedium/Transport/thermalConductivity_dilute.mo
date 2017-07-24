within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_dilute
  "Return thermal conductivity dilute (i.e. density-independent) contribution"
  input ThermodynamicState state;
  output ThermalConductivity lambda_0;

protected
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
  constant Real kilo = 1e3;
  EoS.HelmholtzDerivs f;
  Real cint;
  Real eta_0; // dilute contribution only

algorithm
  // Modelica.Utilities.Streams.print("thermalConductivity_dilute: d = " + String(state.d) + " and T = " + String(state.T));
  tau := state.T/T_red_0;

  // numerator terms, check last numerator exponent
  if (A_num[nDilute_num,2]>-90) then
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num);
  else
    // flagged algorithm, exclude last term from sum loop
    lambda_0 := sum(A_num[i, 1]*tau^A_num[i, 2] for i in 1:nDilute_num-1);
    if (abs(A_num[nDilute_num,2]+99)<eps) then
      // Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 99");
      // used by CO2
      // remember: RefProp uses molar units and g/mol, divide by MolarMass
      f := EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);
      cint := EoS.cp0(f)*f.MM -2.5*f.R*f.MM;
      cint := 1.0 + A_num[nDilute_num, 1]*cint;
      lambda_0 := lambda_0*cint;
    elseif (abs(A_num[nDilute_num,2]+98)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 98");
      assert(false, "not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+97)<eps) then
      Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 97");
      assert(false, "not yet implemented");
    elseif (abs(A_num[nDilute_num,2]+96)<eps) then
      // Modelica.Utilities.Streams.print("thermalConductivity_dilute: flag 96");
      // used by Pentane, Isopentane
      // remember: RefProp uses molar units and g/mol!
      f := EoS.setHelmholtzDerivsSecond(d=state.d,T=state.T);
      cint := EoS.cp0(f) /f.R-2.5;
      eta_0 := dynamicViscosity_dilute(state);
      lambda_0 := (lambda_0*cint+15.0/4.0)*f.R*eta_0/kilo;
    end if;
  end if;

  if (nDilute_den>0) then
    // Modelica.Utilities.Streams.print("thermalConductivity_dilute: use denominator term");
    // used by CO2
    denom := sum(A_den[i, 1]*tau^A_den[i, 2] for i in 1:nDilute_den);
  else
    denom := 1;
  end if;

  lambda_0 := lambda_0/denom*lambda_red_0;
  // Modelica.Utilities.Streams.print(" lambda_0 = " + String(lambda_0));

end thermalConductivity_dilute;
