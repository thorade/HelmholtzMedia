within HelmholtzMedia.HelmholtzFluids;
package R32_RES "R32 with viscosity and thermal conductivity from residual entropy scaling (for comparison)"
  extends HelmholtzMedia.HelmholtzFluids.R32(
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsR32_RES,
    thermalConductivityCoefficients=thermalConductivityCoefficientsR32_RES);

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsR32_RES(
    dynamicViscosityModel=DynamicViscosityModel.RES,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1 "not used",
    epsilon_kappa=1 "not used",
    RES_dilute={1.172179e-10, -1.84014e-07, 9.194887e-05, 0.02459583, 1.0596}
      "dilute gas viscosity in micro Pa s, fitted from 136.34 K to 435.0 K",
    RES_n={0.015284, 0.220607, -0.095156} "fluid-specific RES parameters",
    RES_xi=1) "Viscosity from residual entropy scaling (Martinek et al. 2025)";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsR32_RES(
    thermalConductivityModel=ThermalConductivityModel.RES,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingMolarVolume_background=1 "not used",
    RES_dilute={4.675106e-27, -4.668571e-24, 2.06107e-07, -5.528007e-05, 0.0106548}
      "dilute gas thermal conductivity in W/(m K), fitted from 136.34 K to 435.0 K",
    RES_n={0.83301, 0.496007, 1.424641, -0.174255} "fluid-specific RES parameters",
    RES_xi=1,
    R0=1.03,
    xi_0=1.94e-10,
    Gamma_0=0.0496,
    qd_inverse=5.5829e-10,
    T_ref=526.883) "Thermal conductivity from residual entropy scaling (Li et al. 2024), critical enhancement parameters from the same source";

  annotation (Documentation(info="<html>
<p>
This package is identical to <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.R32\">R32</a>
(equation of state of Tillner-Roth and Yokozeki (1997), the REFPROP 10 default), except for the transport properties:
viscosity and thermal conductivity are computed by residual entropy scaling (RES) with the fluid-specific
parameters of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity, including the parameters of the
critical enhancement). It exists to compare RES with the reference correlations of R32, see
<a href=\"modelica://HelmholtzMedia.Validation.Properties.TransportRESComparison\">TransportRESComparison</a>.
For applications, use <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.R32\">R32</a>.
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
end R32_RES;
