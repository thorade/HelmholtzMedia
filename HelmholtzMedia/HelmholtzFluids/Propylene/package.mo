within HelmholtzMedia.HelmholtzFluids;
package Propylene "Propylene (propene, R-1270) with IIR reference state"
  extends Interfaces.PartialHelmholtzMedium(
    mediumName="propylene" "short name",
    fluidConstants={fluidConstantsPropylene},
    helmholtzCoefficients=helmholtzCoefficientsPropylene,
    thermalConductivityCoefficients=thermalConductivityCoefficientsPropylene,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsPropylene,
    surfaceTensionCoefficients=surfaceTensionCoefficientsPropylene,
    ancillaryCoefficients=ancillaryCoefficientsPropylene,
    fluidLimits=fluidLimitsPropylene,
    Density(min=fluidLimitsPropylene.DMIN, max=fluidLimitsPropylene.DMAX, start=fluidConstantsPropylene.molarMass/fluidConstantsPropylene.criticalMolarVolume),
    Temperature(min=fluidLimitsPropylene.TMIN, max=fluidLimitsPropylene.TMAX, start=298.15),
    AbsolutePressure(min=0, max=fluidLimitsPropylene.PMAX, start=101325),
    SpecificEnthalpy(min=fluidLimitsPropylene.HMIN, max=fluidLimitsPropylene.HMAX, start=(fluidLimitsPropylene.HMIN+fluidLimitsPropylene.HMAX)/2),
    SpecificEntropy(min=fluidLimitsPropylene.SMIN, max=fluidLimitsPropylene.SMAX, start=(fluidLimitsPropylene.SMIN+fluidLimitsPropylene.SMAX)/2));

  final constant FluidConstants
  fluidConstantsPropylene(
    chemicalFormula="C3H6",
    structureFormula="CH2=CH-CH3",
    casRegistryNumber="115-07-1",
    iupacName="propene",
    molarMass=0.04207974,
    hasCriticalData=true,
       criticalTemperature=364.211,
       criticalPressure=4554987.6609,
       criticalMolarVolume=0.04207974/229.62914118 "molar mass/critical density",
       HCRIT0=531088.842,
       SCRIT0=1985.73678,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.366,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=false,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.146,
    triplePointTemperature=87.953,
    triplePointPressure=7.471e-4,
    normalBoilingPoint=225.531,
    meltingPoint=87.953) "Fluid Constants";

  final constant EoSLimits
  fluidLimitsPropylene(
    TMIN=fluidConstantsPropylene.triplePointTemperature,
    TMAX=575,
    DMIN=Modelica.Constants.eps,
    DMAX=23.1*42.07974 "23.1 mol/l",
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-250e3,
    HMAX=1500e3,
    SMIN=-1500,
    SMAX=7000) "Helmholtz EoS Limits";

  final constant HelmholtzMedia.Interfaces.Types.HelmholtzCoefficients
  helmholtzCoefficientsPropylene(
    T_reducing=364.211,
    d_reducing=5.457*42.07974 "5.457 mol/l",
    idealLog=[
       3.0,                 1.0],
    idealPower=[
      -5.18233367588693,    0.0;
       4.36399341249517,    1.0] "IIR reference state",
    idealEinstein=[
       1.544,     -324.0/364.211;
       4.013,     -973.0/364.211;
       8.923,    -1932.0/364.211;
       6.02,     -4317.0/364.211],
    residualPoly=[
       0.04341002,    1.0,     4;
       1.136592,      0.205,   1;
      -0.8528611,     0.56,    1;
       0.5216669,     0.676,   2;
      -1.382953,      1.0,     2;
       0.1214347,     0.5,     3],
    residualBwr=[
      -0.5984662,     1.0,     1,   1;
      -1.391883,      1.94,    1,   2;
      -1.008434,      2.0,     3,   2;
       0.1961249,     1.0,     2,   1;
      -0.360693,      2.66,    2,   2;
      -0.002407175,   0.83,    8,   1],
    residualGauss=[
       0.7432121,     1.6,     1,   2,   2,    -1.07,     -0.77,    1.21,   0.78;
       0.1475162,     2.5,     1,   2,   2,    -0.66,     -0.83,    1.08,   0.82;
      -0.02503391,    3.0,     2,   2,   2,    -1.2,      -0.607,   0.83,   1.94;
      -0.2734409,     2.5,     3,   2,   2,    -1.12,     -0.4,     0.56,   0.69;
       0.006378889,   2.72,    3,   2,   2,    -1.47,     -0.66,    1.22,   1.96;
       0.0150294,     4.0,     2,   2,   2,    -1.93,     -0.07,    1.81,   1.3;
      -0.03162971,    4.0,     1,   2,   2,    -3.3,      -3.1,     1.54,   0.38;
      -0.04107194,    1.0,     2,   2,   2,   -15.4,    -387.0,     1.12,   0.91;
      -1.190241,      4.0,     3,   2,   2,    -6.0,     -41.0,     1.4,    0.7] "<n, t, d, 2, 2, eta, beta, gamma, epsilon>")
  "Coefficients of the Helmholtz EoS";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsPropylene(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=364.211,
    reducingThermalConductivity_0=0.001 "dilute gas correlation in mW/(m K)",
    lambda_0_num_coeffs=[
      -1.37218,     0;
      17.3386,      1;
      -3.27682,     2;
       9.34452,     3;
      12.88,        4;
      -1.5705,      5] "numerator <coeff, power of T/T_red>",
    lambda_0_den_coeffs=[
       1.393679,    0;
      -1.04648,     1;
       1.0,         2] "denominator <coeff, power of T/T_red>",
    reducingTemperature_background=364.211,
    reducingMolarVolume_background=1/5457 "1/(5.457 mol/l)",
    reducingThermalConductivity_background=1.0,
    lambda_b_coeffs=[
      0.00271511,   0,   1,   0;
     -0.0363839,    0,   2,   0;
      0.106159,     0,   3,   0;
     -0.0616755,    0,   4,   0;
      0.0105424,    0,   5,   0;
      0.00994697,   1,   1,   0;
      0.0242705,    1,   2,   0;
     -0.0659429,    1,   3,   0;
      0.0379916,    1,   4,   0;
     -0.0056912,    1,   5,   0] "<coeff, power of T/T_red, power of d/d_red, 0>",
    R0=1.02,
    z=0.065,
    xi_0=0.198e-9,
    Gamma_0=0.057,
    qd_inverse=0.43e-9,
    T_ref=546.32) "Coefficients for the thermal conductivity";

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsPropylene(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=0.4678 "Lennard-Jones size parameter [nm] (from the ECS model of REFPROP 10, not used)",
    epsilon_kappa=298.9 "Lennard-Jones energy parameter [K] (from the ECS model of REFPROP 10, not used)",
    CET=fill(0.0, 0, 2) "no viscosity correlation available")
  "No viscosity correlation available (REFPROP 10 uses extended corresponding states)";

  final constant HelmholtzMedia.Interfaces.Types.SurfaceTensionCoefficients
  surfaceTensionCoefficientsPropylene(
    coeffs=[
      0.05268,    1.186]) "Coefficients for the surface tension <sigma0 and n>";

  final constant HelmholtzMedia.Interfaces.Types.AncillaryCoefficients
  ancillaryCoefficientsPropylene(
    T_crit=364.211,
    p_crit=4555e3,
    d_crit=5.457*42.07974,
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -6.75625,    1.0;
       2.027,      1.5;
      -1.35883,    1.9;
      -2.74671,    4.3;
      -0.936445,  15.0],
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       0.405430,   0.195;
       2.02481,    0.47;
       0.304022,   2.25;
       0.179159,   8.0],
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -1.59841,    0.309;
      -4.73840,    0.853;
     -10.8886,     2.37;
     -31.0312,     5.2;
     -56.9431,    10.0;
    -143.544,     20.0],
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=87.953,
    p_reducing=7.471e-4,
    pressureMelting1=[
      -6593000000.0,   0.0;
       6593000001.0,   2.821],
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL1, DV3, ML1)";

  annotation (Documentation(info="<html>
<p>
These are the coefficients for propylene (propene, R-1270).
The reference state is set to IIR (saturated liquid at 0 &deg;C: h = 200 kJ/kg, s = 1 kJ/(kg K)), the REFPROP default.
</p>
<p>
<b>No viscosity correlation is available:</b> REFPROP 10 computes the viscosity of propylene with an extended corresponding
states model (propane reference fluid), which is not implemented in HelmholtzMedia; <code>dynamicViscosity</code> raises an error.
Since the critical enhancement of the thermal conductivity requires the viscosity, <code>thermalConductivity</code>
returns the dilute gas and background contributions only: compared with REFPROP 10, it is too low by tenths of a percent
in the liquid and gas, by 5.7&nbsp;% at 350 K and 385 kg/m3 and by up to tens of percent near the critical point
(36&nbsp;% at 370 K and 230 kg/m3). The transport properties of propylene will be completed later.
</p>
<p>
The ancillary equations (vapour pressure, saturated densities) are those of Lemmon et al. as given in REFPROP 10,
the melting line is that of Reeves et al. (1964) with the triple point pressure of Lemmon et al.
</p>

<dl>
<dt> Lemmon, E.W., Overhoff, U., McLinden, M.O. and Wagner, W.</dt>
<dd> <b>A Reference Equation of State for the Thermodynamic Properties of Propene for Temperatures from the Melting Line to 575 K and Pressures up to 1000 MPa</b><br />
     to be submitted to J. Phys. Chem. Ref. Data; coefficients as implemented in REFPROP 10.
</dd>
<dt> Assael, M.J., Koutian, A., Huber, M.L. and Perkins, R.A.</dt>
<dd> <b>Reference Correlations of the Thermal Conductivity of Ethene and Propene</b><br />
     J. Phys. Chem. Ref. Data, 45:033104, 2016.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4958984\">10.1063/1.4958984</a>
</dd>
<dt> Mulero, A., Cachadi&ntilde;a, I. and Parra, M.I.</dt>
<dd> <b>Recommended Correlations for the Surface Tension of Common Fluids</b><br />
     J. Phys. Chem. Ref. Data, 41:043105, 2012.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4768782\">10.1063/1.4768782</a>
</dd>
<dt> Reeves, L.E., Scott, G.J. and Babb, S.E.</dt>
<dd> <b>Melting Curves of Pressure-Transmitting Fluids</b><br />
     J. Chem. Phys., 40:3662-3666, 1964.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.1725068\">10.1063/1.1725068</a>
</dd>
<dt> Lemmon, E.W., Bell, I.H., Huber, M.L. and McLinden, M.O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP, Version 10.0</b><br />
     National Institute of Standards and Technology, 2018.<br />
     DOI: <a href=\"https://doi.org/10.18434/T4/1502528\">10.18434/T4/1502528</a>
</dd>
</dl>
</html>"));
end Propylene;
