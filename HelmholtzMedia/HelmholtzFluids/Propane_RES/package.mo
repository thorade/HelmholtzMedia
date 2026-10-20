within HelmholtzMedia.HelmholtzFluids;
package Propane_RES "Propane with viscosity and thermal conductivity from residual entropy scaling (for comparison)"
  extends HelmholtzMedia.HelmholtzFluids.Propane(
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsPropane_RES,
    thermalConductivityCoefficients=thermalConductivityCoefficientsPropane_RES);

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsPropane_RES(
    dynamicViscosityModel=DynamicViscosityModel.RES,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1 "not used",
    epsilon_kappa=1 "not used",
    RES_dilute={5.425586e-12, -1.475908e-08, 5.310593e-06, 0.02684617, 8.35646e-15}
      "dilute gas viscosity in micro Pa s, fitted from 85.525 K to 650.0 K",
    RES_n={0.406055, -0.178357, 0.041571} "fluid-specific RES parameters",
    RES_xi=1) "Viscosity from residual entropy scaling (Martinek et al. 2025)";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsPropane_RES(
    thermalConductivityModel=ThermalConductivityModel.RES,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingMolarVolume_background=1 "not used",
    RES_dilute={2.784026e-27, -4.412509e-24, 1.457529e-07, 2.207303e-05, -0.00124778}
      "dilute gas thermal conductivity in W/(m K), fitted from 85.525 K to 650.0 K",
    RES_n={8.726474, -13.787569, 8.781108, -1.303329} "fluid-specific RES parameters",
    RES_xi=1,
    R0=1.03,
    xi_0=1.94e-10,
    Gamma_0=0.0496,
    qd_inverse=7.1664e-10,
    T_ref=554.73) "Thermal conductivity from residual entropy scaling (Li et al. 2024), critical enhancement parameters from the same source";

  annotation (Documentation(info="<html>
<p>
This package is identical to <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.Propane\">Propane</a>
(equation of state of Lemmon et al. (2009), the REFPROP 10 default), except for the transport properties:
viscosity and thermal conductivity are computed by residual entropy scaling (RES) with the fluid-specific
parameters of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity, including the parameters of the
critical enhancement). It exists to compare RES with the reference correlations of Propane, see
<a href=\"modelica://HelmholtzMedia.Validation.Properties.TransportRESComparison\">TransportRESComparison</a>.
For applications, use <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.Propane\">Propane</a>.
</p>
<p>
The reference state is set to IIR (saturated liquid at 0 &deg;C: h = 200 kJ/kg, s = 1 kJ/(kg K)), the REFPROP default.
</p>
<dl>
<dt> Martinek, V., Bell, I.H., Herzog, R., Richter, M. and Yang, X.</dt>
<dd> <b>Entropy Scaling of Viscosity IV - Application to 124 Industrially Important Fluids</b><br />
     J. Chem. Eng. Data, 70:727-742, 2025.<br />
     DOI: <a href=\"https://doi.org/10.1021/acs.jced.4c00451\">10.1021/acs.jced.4c00451</a>
</dd>
<dt> Li, Z., Duan, Y. and Yang, X.</dt>
<dd> <b>Linking Thermal Conductivity to Equations of State Using the Residual Entropy Scaling Theory</b><br />
     Ind. Eng. Chem. Res., 63(42):18160-18175, 2024.<br />
     DOI: <a href=\"https://doi.org/10.1021/acs.iecr.4c02946\">10.1021/acs.iecr.4c02946</a>
</dd>
</dl>
</html>"));
end Propane_RES;
