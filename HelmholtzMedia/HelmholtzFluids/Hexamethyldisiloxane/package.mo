within HelmholtzMedia.HelmholtzFluids;
package Hexamethyldisiloxane "hexamethyldisiloxane (MM)"
extends Interfaces.PartialHelmholtzMedium(
  fluidConstants={fluidConstantsMM},
  helmholtzCoefficients=helmholtzCoefficientsMM,
  thermalConductivityCoefficients=thermalConductivityCoefficientsMM,
  dynamicViscosityCoefficients=dynamicViscosityCoefficientsMM,
  surfaceTensionCoefficients=surfaceTensionCoefficientsMM,
  ancillaryCoefficients=ancillaryCoefficientsMM,
  fluidLimits=fluidLimitsMM,
  Density(min=fluidLimitsMM.DMIN, max=fluidLimitsMM.DMAX, start=fluidConstantsMM.molarMass/fluidConstantsMM.criticalMolarVolume),
  Temperature(min=fluidLimitsMM.TMIN, max=fluidLimitsMM.TMAX, start=298.15),
  AbsolutePressure(min=0, max=fluidLimitsMM.PMAX, start=101325),
  SpecificEnthalpy(min=fluidLimitsMM.HMIN, max=fluidLimitsMM.HMAX, start=0),
  SpecificEntropy(min=fluidLimitsMM.SMIN, max=fluidLimitsMM.SMAX, start=0));

  final constant FluidConstants
  fluidConstantsMM(
    chemicalFormula="C6H18OSi2",
    structureFormula="",
    casRegistryNumber="107-46-0",
    iupacName="",
    molarMass=0.1623768,
    hasCriticalData=true,
       criticalTemperature=518.7,
       criticalPressure=1931130.604 "from the EoS at T_crit and d_crit",
       criticalMolarVolume=1/1653,
       HCRIT0=376537.80456143,
       SCRIT0=829.07069762325,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.801,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.418,
    triplePointTemperature=204.93,
    triplePointPressure=2.955 "from the EoS at the triple point temperature",
    normalBoilingPoint=373.658,
    meltingPoint=204.93) "Fluid Constants";

  final constant EoSLimits
  fluidLimitsMM(
    TMIN=fluidConstantsMM.triplePointTemperature,
    TMAX=673 "range of validity of the EoS up to 580 K, extrapolated above",
    DMIN=Modelica.Constants.eps,
    DMAX=930,
    PMIN=Modelica.Constants.small,
    PMAX=130e6,
    HMIN=-320e3,
    HMAX=835e3,
    SMIN=-1260,
    SMAX=3630) "Helmholtz EoS Limits";

  final constant HelmholtzMedia.Interfaces.Types.HelmholtzCoefficients
  helmholtzCoefficientsMM(
    T_reducing=518.7,
    d_reducing=0.1623768*1653 "1.653 mol/l",
    idealLog=[
       3.0,                 1.0],
    idealPower=[
      72.110753778474722,   0.0;
     -10.43149917525645,    1.0] "NBP reference state",
    idealEinstein=[
      18.59,       -20.0/518.7;
      29.58,     -1400.0/518.7;
      19.74,     -3600.0/518.7;
       4.87,     -6300.0/518.7],
    residualPoly=[
       0.05063651,    1.0,     4;
       8.604724,      0.346,   1;
      -9.179684,      0.46,    1;
      -1.146325,      1.01,    2;
       0.4878559,     0.59,    3],
    residualBwr=[
      -2.434088,      2.6,     1,   2;
      -1.621326,      3.33,    3,   2;
       0.6239872,     0.75,    2,   1;
      -2.306057,      2.95,    2,   2;
      -0.05555096,    0.93,    7,   1],
    residualGauss=[
       9.385015,      1.33,    1,   2,   2,    -1.0334,   -0.4707,   1.7754,   0.8927;
      -2.493508,      1.68,    1,   2,   2,    -1.544,    -0.32,     0.692,    0.5957;
      -3.308032,      1.7,     3,   2,   2,    -1.113,    -0.404,    1.242,    0.559;
      -0.1885803,     3.08,    3,   2,   2,    -1.113,    -0.517,    0.421,    1.056;
      -0.09883865,    5.41,    1,   2,   2,    -1.11,     -0.432,    0.406,    1.3;
       0.111109,      1.4,     2,   2,   2,    -7.2,      -7.2,      0.163,    0.106;
       0.1061928,     1.1,     3,   2,   2,    -1.45,     -1.2,      0.795,    0.181;
      -0.01452454,    5.3,     1,   2,   2,    -4.73,    -35.8,      0.88,     0.525] "<n, t, d, 2, 2, eta, beta, gamma, epsilon>")
  "Coefficients of the Helmholtz EoS";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsMM(
    thermalConductivityModel=ThermalConductivityModel.RES,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingMolarVolume_background=1,
    RES_dilute={2.069178e-13, -4.626304e-10, 3.990619e-07, -6.735081e-05, 7.961889e-03}
      "dilute gas thermal conductivity in W/(m K), fitted from 204.93 K to 580 K",
    RES_n={0, -2.998729, 4.897334, -0.700838} "fluid-specific RES parameters",
    RES_xi=1,
    R0=1.02,
    xi_0=2.68e-10,
    Gamma_0=0.062,
    qd_inverse=8.4e-10,
    T_ref=778.05) "Thermal conductivity from residual entropy scaling (Li et al. 2024), critical enhancement parameters from the same source";

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsMM(
    dynamicViscosityModel=DynamicViscosityModel.RES,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=1,
    epsilon_kappa=1,
    RES_dilute={3.045015e-11, -5.987369e-08, 3.713322e-05, 1.303627e-02, 6.489577e-01}
      "dilute gas viscosity in micro Pa s, fitted from 204.93 K to 580 K",
    RES_n={0.575147, -0.299109, 0.069042} "fluid-specific RES parameters",
    RES_xi=1) "Viscosity from residual entropy scaling (Martinek et al. 2025)";

  final constant HelmholtzMedia.Interfaces.Types.SurfaceTensionCoefficients
  surfaceTensionCoefficientsMM(
    coeffs=[
      0.04576,    1.272]) "Coefficients for the surface tension (Mulero and Cachadina 2014)";

  final constant HelmholtzMedia.Interfaces.Types.AncillaryCoefficients
  ancillaryCoefficientsMM(
    T_crit=518.7,
    p_crit=1931130,
    d_crit=0.1623768*1653,
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -8.5023,     1.0;
       3.803,      1.5;
      -3.415,      1.98;
      -4.679,      3.86;
      -3.106,     14.6],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       4.003,      0.436;
      -6.406,      0.827;
      11.5,        1.24;
     -10.04,       1.7;
       4.0,        2.23],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -3.7421,     0.428;
     -37.087,      1.79;
      75.46,       2.28;
     -71.67,       2.8;
     -68.69,       7.0;
    -178.4,       15.4])
  "Coefficients for the ancillary equations (PS5, DL1, DV3) of Thol et al. (2016)";


  annotation (Documentation(info="<html>
These are the coefficients for MM.
There is no dedicated viscosity or thermal conductivity correlation for MM (REFPROP 10 uses an extended corresponding
states model with nitrogen as reference fluid). Both are computed by residual entropy scaling (RES) with the fluid-specific
parameters of Martinek et al. (2025, viscosity) and Li et al. (2024, thermal conductivity, including the parameters of the critical enhancement).
The RES parameters were fitted with the residual entropy of the same equation of state (Thol et al. 2016, the REFPROP 10 default).
The equation of state is that of Thol et al. (2016), the REFPROP 10 default, with the ancillary equations given there
(valid from the triple point at 204.93 K to 580 K and up to 130 MPa; the temperature limit TMAX = 673 K extrapolates).
The reference state is set to NBP (saturated liquid at 101325 Pa: h = 0, s = 0), the REFPROP default.

<dl>
<dt> Thol, M., Dubberke, F.H., Rutkai, G., Windmann, T., K&ouml;ster, A., Span, R. and Vrabec, J.</dt>
<dd> <b>Fundamental Equation of State Correlation for Hexamethyldisiloxane Based on Experimental and Molecular Simulation Data</b><br />
     Fluid Phase Equilibria, 418:133-151, 2016<br />
     DOI: <a href=\"https://doi.org/10.1016/j.fluid.2015.09.047\">10.1016/j.fluid.2015.09.047</a>
</dd>
<dt> Martinek, V., Bell, I.H., Herzog, R., Richter, M. and Yang, X.</dt>
<dd> <b>Entropy Scaling of Viscosity IV - Application to 124 Industrially Important Fluids</b><br />
     J. Chem. Eng. Data, 2025.
</dd>
<dt> Li, Z., Duan, Y. and Yang, X.</dt>
<dd> <b>Linking Thermal Conductivity to Equations of State Using the Residual Entropy Scaling Theory</b><br />
     Ind. Eng. Chem. Res., 63(42):18160-18175, 2024.
</dd>
<dt> Mulero, A. and Cachadi&ntilde;a, I.</dt>
<dd> <b>Recommended Correlations for the Surface Tension of Several Fluids Included in the REFPROP Program</b><br />
     J. Phys. Chem. Ref. Data, 43:023104, 2014.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4878755\">10.1063/1.4878755</a>
</dd>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
</dl>
</html>"));

end Hexamethyldisiloxane;
