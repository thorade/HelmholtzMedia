within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_background
  "Return thermal conductivity background/residual  contribution"
  input ThermodynamicState state;
  output ThermalConductivity lambda_b;

protected
  ThermalConductivityModel thermalConductivityModel=thermalConductivityCoefficients.thermalConductivityModel;
  MolarMass MM = fluidConstants[1].molarMass;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_background=fluidConstants[1].molarMass/thermalConductivityCoefficients.reducingMolarVolume_background;
  Real delta "reduced density";
  Temperature T_red_background=thermalConductivityCoefficients.reducingTemperature_background;
  Real tau "reduced temperature";

  // coeffs for background contribution
  Integer nBackground = size(thermalConductivityCoefficients.lambda_b_coeffs, 1);
  Real[nBackground, 4] B= thermalConductivityCoefficients.lambda_b_coeffs;

  Real lambda_red_background=thermalConductivityCoefficients.reducingThermalConductivity_background;

  // helium (Hands and Arp 1981), only used for TC0
  Temperature T;
  Density d;
  Real lambda_0;
  Real lambda_e;
  Real lambda_c;
  constant Real c[11]={1.862970530e-4, -7.275964435e-7, -1.427549651e-4, 3.290833592e-5, -5.213335363e-8, 4.492659933e-8,
                       -5.924416513e-9, 7.087321137e-6, -6.013335678e-6, 8.067145814e-7, 3.995125013e-7};
  constant Real x0=0.392;
  constant Real E1=2.8461;
  constant Real E2=0.27156;
  constant Real beta=0.3554;
  constant Real gamma=1.1743;
  constant Real delta_He=4.304;
  constant Real rhoc_He=69.158;
  constant Real Tc_He=5.18992;
  constant Real pc_He=2.2746e5;
  Real DeltaT;
  Real DeltaRho;
  EoS.HelmholtzDerivs f;
  Real dpdT_d;
  Real K_T;
  Real K_Tprime;
  Real K_Tbar;
  Real W;
  Real x;
  Real y;
  Real h;
  Real dhdx;

algorithm
  // initialize output
  lambda_b := 0;

  if (thermalConductivityModel == ThermalConductivityModel.TC0) then
    // hardcoded models, use mediumName to distinguish further
    if mediumName == "helium" then
      // Hands and Arp (1981) as implemented in REFPROP 10 (TC7 in HELIUM.FLD); d in kg/m3, lambda in W/(m K)
      T := state.T;
      d := state.d;
      lambda_0 := 2.7870034e-3*T^0.7034007057*exp(3.739232544/T - 2.620316969e1/T^2 + 5.982252246e1/T^3 - 4.926397634e1/T^4);
      lambda_e := (c[1] + c[2]*T + c[3]*T^(1/3) + c[4]*T^(2/3))*d
                + (c[5] + c[6]*T^(1/3) + c[7]*T^(2/3))*d^3
                + (c[8] + c[9]*T^(1/3) + c[10]*T^(2/3) + c[11]/T)*d^2*log(d/68.0);
      lambda_c := 0;
      if (T > 3.5) and (T < 12) then
        DeltaT := abs(1 - T/Tc_He);
        DeltaRho := max(abs(1 - d/rhoc_He), 1e-12);
        f := EoS.setHelmholtzDerivsSecond(T=T, d=d, phase=1);
        dpdT_d := EoS.dpTd(f);
        K_T := 1/(d*EoS.dpdT(f));
        W := (DeltaT/0.2)^2 + (DeltaRho/0.25)^2;
        if (W > 1) then
          K_Tbar := K_T;
        else
          // scaling variable x = DeltaT/DeltaRho^(1/beta), as in the REFPROP 10 implementation
          x := DeltaT/DeltaRho^(1/beta);
          y := 1 + E2*(1 + x/x0)^(2*beta);
          h := E1*(1 + x/x0)*y^((gamma - 1)/(2*beta));
          dhdx := E1/x0*y^((gamma - 1)/(2*beta))*(1 + (gamma - 1)*E2*(1 + x/x0)^(2*beta)/y);
          K_Tprime := 1/(DeltaRho^(delta_He - 1)*(delta_He*h - x/beta*dhdx)*(d/rhoc_He)^2*pc_He);
          K_Tbar := W*K_T + (1 - W)*K_Tprime;
        end if;
        // the factors 3.4685233e-17 and 3.726229668 are taken from REFPROP (converted to SI units)
        lambda_c := 3.4685233e-17*3.726229668*sqrt(K_Tbar)*T^2/d/dynamicViscosity(state)*dpdT_d^2
                    *exp(-18.66*DeltaT^2 - 4.25*DeltaRho^4);
      end if;
      lambda_b := lambda_0 + lambda_e + lambda_c;
    else
      assert(false, "thermalConductivity: no hardcoded TC0 model for " + mediumName);
    end if;
  elseif (thermalConductivityModel == ThermalConductivityModel.TC1) then
    tau := state.T/T_red_background;
    delta := state.d/d_red_background;
    lambda_b := sum((B[i, 1]*tau^B[i, 2])*(delta)^B[i, 3] for i in 1:nBackground);
    lambda_b := lambda_b*lambda_red_background;
  elseif (thermalConductivityModel == ThermalConductivityModel.TC2) then
    assert(false, "ThermalconductivityModel TC2 not yet implemented");
  else
    assert(false, "unknown ThermalconductivityModel");
  end if;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print(" lambda_b = " + String(lambda_b));
  Modelica.Utilities.Streams.print("===========================================");
  */

  annotation (Documentation(info="<html>
  <p>
The thermal conductivity (TC) is split into three parts: ideal gas TC lamda_0, residual TC lambda_r and critical TC enhancement lambda_c.
Sometimes the residual TC is again split into two parts.
This allows to develop functions for each contribution separately.
The sum of ideal gas TC and residual TC is called background TC.
Ideal gas TC depends on Temperature only and can be modelled by a quadratic function.
Residual TC is also modeled by a polynominal.
At the critical point TC becomes infinite; TC is enhanced for a large region around the critical point.
The critical enhancement can be described by various alternative approaches.
Here, the simplified approach as suggested by Olchowy and Sengers is implemented.

Special thanks go to Eric W. Lemmon for answering all my emails
and programming a special version of RefProp that outputs also intermediate values.

</p>
<dl>
<dt>Olchowy, G.A. and Sengers, J.V</dt>
<dd> <b>A simplified representation for the thermal conductivity of fluids in the critical region</b>.<br />
     International Journal of Thermophysics (1998) 10, 417-426.<br />
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end thermalConductivity_background;
