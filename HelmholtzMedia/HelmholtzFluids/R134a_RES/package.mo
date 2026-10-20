within HelmholtzMedia.HelmholtzFluids;
package R134a_RES "R134a with viscosity and thermal conductivity from residual entropy scaling (for comparison)"
  extends HelmholtzMedia.HelmholtzFluids.R134a(
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsR134a_RES,
    thermalConductivityCoefficients=thermalConductivityCoefficientsR134a_RES);

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsR134a_RES(
    dynamicViscosityModel=DynamicViscosityModel.RES,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1 "not used",
    epsilon_kappa=1 "not used",
    RES_dilute={-1.269154e-12, 4.4985e-09, -9.057271e-06, 0.04267365, -0.1766998}
      "dilute gas viscosity in micro Pa s, fitted from 169.85 K to 455.0 K",
    RES_n={0.310059, -0.106155, 0.020416} "fluid-specific RES parameters",
    RES_xi=1) "Viscosity from residual entropy scaling (Martinek et al. 2025)";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsR134a_RES(
    thermalConductivityModel=ThermalConductivityModel.RES,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingMolarVolume_background=1 "not used",
    RES_dilute={-2.822913e-26, 3.607363e-23, -1.701545e-20, 8.00982e-05, -0.0105248}
      "dilute gas thermal conductivity in W/(m K), fitted from 169.85 K to 455.0 K",
    RES_n={1.388264, -0.760952, 1.664978, -0.125953} "fluid-specific RES parameters",
    RES_xi=1,
    R0=1.03,
    xi_0=1.94e-10,
    Gamma_0=0.0496,
    qd_inverse=5.2854e-10,
    T_ref=561.411) "Thermal conductivity from residual entropy scaling (Li et al. 2024), critical enhancement parameters from the same source";

  annotation (Documentation(info="<html>
<p>
This package is identical to <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.R134a\">R134a</a>
(equation of state of Tillner-Roth and Baehr (1994), the REFPROP 10 default), except for the transport properties:
viscosity and thermal conductivity are computed by residual entropy scaling (RES) with the fluid-specific
parameters of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity, including the parameters of the
critical enhancement). It exists to compare RES with the reference correlations of R134a, see
<a href=\"modelica://HelmholtzMedia.Validation.Properties.TransportRESComparison\">TransportRESComparison</a>.
For applications, use <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.R134a\">R134a</a>.
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
end R134a_RES;
