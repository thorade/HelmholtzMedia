within HelmholtzMedia.HelmholtzFluids;
package Carbondioxide_RES "Carbondioxide with viscosity and thermal conductivity from residual entropy scaling (for comparison)"
  extends HelmholtzMedia.HelmholtzFluids.Carbondioxide(
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsCarbondioxide_RES,
    thermalConductivityCoefficients=thermalConductivityCoefficientsCarbondioxide_RES);

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsCarbondioxide_RES(
    dynamicViscosityModel=DynamicViscosityModel.RES,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1 "not used",
    epsilon_kappa=1 "not used",
    RES_dilute={-1.425519e-12, 9.864163e-09, -2.85092e-05, 0.06274267, -1.484764}
      "dilute gas viscosity in micro Pa s, fitted from 216.592 K to 2000.0 K",
    RES_n={0.2335, 0.007837, -0.023471} "fluid-specific RES parameters",
    RES_xi=1) "Viscosity from residual entropy scaling (Martinek et al. 2025)";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsCarbondioxide_RES(
    thermalConductivityModel=ThermalConductivityModel.RES,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingMolarVolume_background=1 "not used",
    RES_dilute={4.004284e-15, -1.789454e-11, 1.587771e-08, 7.555339e-05, -0.006818247}
      "dilute gas thermal conductivity in W/(m K), fitted from 216.592 K to 2000.0 K",
    RES_n={1.706852, -2.225083, 2.920715, -0.376879} "fluid-specific RES parameters",
    RES_xi=1,
    R0=1.02,
    xi_0=1.5e-10,
    Gamma_0=0.052,
    qd_inverse=4.0e-10,
    T_ref=456.19) "Thermal conductivity from residual entropy scaling (Li et al. 2024), critical enhancement parameters from the same source";

  annotation (Documentation(info="<html>
<p>
This package is identical to <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.Carbondioxide\">Carbondioxide</a>
(equation of state of Span and Wagner (1996), the REFPROP 10 default), except for the transport properties:
viscosity and thermal conductivity are computed by residual entropy scaling (RES) with the fluid-specific
parameters of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity, including the parameters of the
critical enhancement). It exists to compare RES with the reference correlations of carbon dioxide, see
<a href=\"modelica://HelmholtzMedia.Validation.Properties.TransportRESComparison\">TransportRESComparison</a>.
For applications, use <a href=\"modelica://HelmholtzMedia.HelmholtzFluids.Carbondioxide\">Carbondioxide</a>.
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
end Carbondioxide_RES;
