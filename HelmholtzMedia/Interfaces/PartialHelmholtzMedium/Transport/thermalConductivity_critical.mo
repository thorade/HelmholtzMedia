within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function thermalConductivity_critical
  "Return thermal conductivity critical enhancement"
  input ThermodynamicState state;
  output ThermalConductivity lambda_c;

protected
  MolarMass MM = fluidConstants[1].molarMass;
  Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  Density d_red_residual=fluidConstants[1].molarMass/
      thermalConductivityCoefficients.reducingMolarVolume_residual;

  Temperature T_crit=fluidConstants[1].criticalTemperature;
  Real tau "reduced temperature";

  AbsolutePressure p_crit=fluidConstants[1].criticalPressure;

  // coeffs for critical enhancement
  Real nu=thermalConductivityCoefficients.nu;
  Real gamma=thermalConductivityCoefficients.gamma;
  Real R0=thermalConductivityCoefficients.R0;
  Real z=thermalConductivityCoefficients.z;
  Real c=thermalConductivityCoefficients.c;
  Real xi_0=thermalConductivityCoefficients.xi_0;
  Real Gamma_0=thermalConductivityCoefficients.Gamma_0;
  Real q_D=1/thermalConductivityCoefficients.qd_inverse;
  Temperature T_ref=thermalConductivityCoefficients.T_ref;

  // interim variables for critical enhancement
  constant Real pi=Modelica.Constants.pi;
  constant Real k_b=Modelica.Constants.k;
  EoS.HelmholtzDerivs f;
  EoS.HelmholtzDerivs f_ref;
  Real ddpT;
  Real ddpT_ref;
  Real chi;
  Real chi_ref;
  Real Delta_chi;
  Real xi;
  Real Omega_0;
  Real Omega;

  SpecificHeatCapacity Cp;
  SpecificHeatCapacity Cv;
  DynamicViscosity eta_b;

algorithm
  // critical enhancement by the simplified crossover model by Olchowy and Sengers
  if ((state.T > T_ref) or (state.d < d_crit/100)) then
    lambda_c := 0; // far away from critical point
  else
    // use critical values from EoS to calculate chi, Omega and lambda_c
    // watch out: algorithm for chi and chi_ref are different (chi_ref is multiplied with T_ref/state.T)
    f     := EoS.setHelmholtzDerivsSecond(T=state.T, d=state.d, phase=1);
    f_ref := EoS.setHelmholtzDerivsSecond(T=T_ref,   d=state.d, phase=1);
    ddpT     := 1.0/EoS.dpdT(f);
    ddpT_ref := 1.0/EoS.dpdT(f_ref);
    chi     := p_crit/d_crit^2*state.d*ddpT;
    chi_ref := p_crit/d_crit^2*state.d*ddpT_ref*T_ref/state.T;

    Delta_chi := chi - chi_ref;

    if (Delta_chi < 0) then
      lambda_c := 0;
    else
      xi := xi_0*(Delta_chi/Gamma_0)^(nu/gamma);

      Cp := specificHeatCapacityCp(state=state);
      Cv := specificHeatCapacityCv(state=state);
      Omega := 2/pi*((Cp - Cv)/Cp*atan(q_D*xi) + Cv/Cp*q_D*xi);
      Omega_0 := 2/pi*(1 - exp(-1/(1/(q_D*xi) + ((q_D*xi*d_crit/state.d)^2)/3)));

      eta_b := dynamicViscosity(state=state);
      lambda_c := (state.d*Cp*R0*k_b*state.T)/(6*pi*eta_b*xi)*(Omega - Omega_0);
      lambda_c := max(0, lambda_c);
    end if;
  end if;

  /* // following lines are for debugging only
  Modelica.Utilities.Streams.print("===========================================");
  Modelica.Utilities.Streams.print("        d = " + String(state.d) + " and T = " + String(state.T));
  Modelica.Utilities.Streams.print("   ddpT   = " + String(ddpT) + " and ddpT_ref = " + String(ddpT_ref));
  Modelica.Utilities.Streams.print("   chi    = " + String(chi) + "  and  chi_ref = " + String(chi_ref));
  Modelica.Utilities.Streams.print("Delta_chi = " + String(Delta_chi));
  Modelica.Utilities.Streams.print("       xi = " + String(xi));
  Modelica.Utilities.Streams.print("       Cp = " + String(Cp) + "  and  Cv = " + String(Cv));
  Modelica.Utilities.Streams.print("  Omega_0 = " + String(Omega_0));
  Modelica.Utilities.Streams.print("    Omega = " + String(Omega));
  Modelica.Utilities.Streams.print("    eta_b = " + String(eta_b));
  Modelica.Utilities.Streams.print(" lambda_c = " + String(lambda_c));
  Modelica.Utilities.Streams.print("===========================================");
  */

  annotation (Documentation(info="<html>
  <p>
The thermal conductivity (TC) is split into three parts: ideal gas TC lamda_0, residual TC lambda_r and critical TC enhancement lambda_c.
Sometimes the residual TC is again split into two parts.
This allows to develop functions for each contribution seperately.
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
<dd> <b>A simplified representation for the thermal conductivity of fluids in the critical region</b>.<br>
     International Journal of Thermophysics (1998) 10, 417-426.<br>
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end thermalConductivity_critical;
