within HelmholtzMedia.HelmholtzFluids;
package Ammonia "Ammonia (R-717) with IIR reference state"
  extends Interfaces.PartialHelmholtzMedium(
    mediumName="ammonia" "short name",
    fluidConstants={fluidConstantsAmmonia},
    helmholtzCoefficients=helmholtzCoefficientsAmmonia,
    thermalConductivityCoefficients=thermalConductivityCoefficientsAmmonia,
    dynamicViscosityCoefficients=dynamicViscosityCoefficientsAmmonia,
    surfaceTensionCoefficients=surfaceTensionCoefficientsAmmonia,
    ancillaryCoefficients=ancillaryCoefficientsAmmonia,
    fluidLimits=fluidLimitsAmmonia,
    Density(min=fluidLimitsAmmonia.DMIN, max=fluidLimitsAmmonia.DMAX, start=fluidConstantsAmmonia.molarMass/fluidConstantsAmmonia.criticalMolarVolume),
    Temperature(min=fluidLimitsAmmonia.TMIN, max=fluidLimitsAmmonia.TMAX, start=298.15),
    AbsolutePressure(min=0, max=fluidLimitsAmmonia.PMAX, start=101325),
    SpecificEnthalpy(min=fluidLimitsAmmonia.HMIN, max=fluidLimitsAmmonia.HMAX, start=(fluidLimitsAmmonia.HMIN+fluidLimitsAmmonia.HMAX)/2),
    SpecificEntropy(min=fluidLimitsAmmonia.SMIN, max=fluidLimitsAmmonia.SMAX, start=(fluidLimitsAmmonia.SMIN+fluidLimitsAmmonia.SMAX)/2));

  final constant FluidConstants
  fluidConstantsAmmonia(
    chemicalFormula="NH3",
    structureFormula="NH3",
    casRegistryNumber="7664-41-7",
    iupacName="ammonia",
    molarMass=0.01703052,
    hasCriticalData=true,
       criticalTemperature=405.56,
       criticalPressure=11363395.0144,
       criticalMolarVolume=0.01703052/233.25000192 "molar mass/critical density",
       HCRIT0=1102127.53,
       SCRIT0=3511.6783,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=1.47,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.256,
    triplePointTemperature=195.49,
    triplePointPressure=6053.39,
    normalBoilingPoint=239.834,
    meltingPoint=195.49) "Fluid Constants";

  final constant EoSLimits
  fluidLimitsAmmonia(
    TMIN=fluidConstantsAmmonia.triplePointTemperature,
    TMAX=725,
    DMIN=Modelica.Constants.eps,
    DMAX=53.13*17.03052 "53.13 mol/l",
    PMIN=Modelica.Constants.small,
    PMAX=1000e6,
    HMIN=-200e3,
    HMAX=3000e3,
    SMIN=-600,
    SMAX=15000) "Helmholtz EoS Limits";

  final constant HelmholtzMedia.Interfaces.Types.HelmholtzCoefficients
  helmholtzCoefficientsAmmonia(
    T_reducing=405.56,
    d_reducing=13.696*17.03052 "13.696 mol/l",
    idealLog=[
       3.0,                 1.0],
    idealPower=[
      -5.60372393904016,    0.0;
       4.86527247654725,    1.0] "IIR reference state",
    idealEinstein=[
       2.224,     -1646.0/405.56;
       3.148,     -3965.0/405.56;
       0.9579,    -7231.0/405.56],
    residualPoly=[
       0.006132232,   1.0,     4;
       1.7395866,     0.382,   1;
      -2.2261792,     1.0,     1;
      -0.30127553,    1.0,     2;
       0.08967023,    0.677,   3],
    residualBwr=[
      -0.076387037,   2.915,   3,   2;
      -0.84063963,    3.51,    2,   2;
      -0.27026327,    1.063,   3,   1],
    residualGauss=[
       6.212578,      0.655,   1,   2,   2,    -0.42776,     -1.708,    1.036,     -0.0726;
      -5.7844357,     1.3,     1,   2,   2,    -0.6424,      -1.4865,   1.2777,    -0.1274;
       2.4817542,     3.1,     1,   2,   2,    -0.8175,      -2.0915,   1.083,      0.7527;
      -2.3739168,     1.4395,  2,   2,   2,    -0.7995,      -2.43,     1.2906,     0.57;
       0.01493697,    1.623,   2,   2,   2,    -0.91,        -0.488,    0.928,      2.2;
      -3.7749264,     0.643,   1,   2,   2,    -0.3574,      -1.1,      0.934,     -0.243;
       0.0006254348,  1.13,    3,   2,   2,    -1.21,        -0.85,     0.919,      2.96;
      -0.000017359,   4.5,     3,   2,   2,    -4.14,        -1.14,     1.852,      3.02;
      -0.13462033,    1.0,     1,   2,   2,   -22.56,      -945.64,    1.05897,    0.9574;
       0.07749072839, 4.0,     1,   2,   2,   -22.68,      -993.85,    1.05277,    0.9576] "<n, t, d, 2, 2, eta, beta, gamma, epsilon>",
    residualAssociating=[
      -1.6909858,     4.3315,  1,    -2.8452,   0.3696,   1.108,   0.4478,   1.244;
       0.93739074,    4.015,   1,    -2.8342,   0.2962,   1.313,   0.44689,  0.6826] "<n, t, d, eta, beta, gamma, epsilon, b>")
  "Coefficients of the Helmholtz EoS";

  final constant HelmholtzMedia.Interfaces.Types.ThermalConductivityCoefficients
  thermalConductivityCoefficientsAmmonia(
    thermalConductivityModel=ThermalConductivityModel.TC1,
    thermalConductivityCriticalEnhancementModel=ThermalConductivityCriticalEnhancementModel.TK3,
    reducingTemperature_0=405.56,
    reducingThermalConductivity_0=0.001 "dilute gas correlation in mW/(m K)",
    lambda_0_num_coeffs=[
       86.9294,     0;
     -170.5502,     1;
      608.0287,     2;
     -100.9764,     3;
       85.1986,     4] "numerator <coeff, power of T/T_red>",
    lambda_0_den_coeffs=[
        4.68994,    0;
        9.21307,    1;
       -1.53637,    2;
        1.0,        3] "denominator <coeff, power of T/T_red>",
    reducingTemperature_background=405.56,
    reducingMolarVolume_background=1/13696 "1/(13.696 mol/l)",
    reducingThermalConductivity_background=1.0,
    lambda_b_coeffs=[
      0.1034320,    0,   1,   0;
     -0.1125970,    0,   2,   0;
      0.2333010,    0,   3,   0;
     -0.1125360,    0,   4,   0;
      0.0141129,    0,   5,   0;
     -0.0283976,    1,   1,   0;
      0.0482520,    1,   2,   0;
     -0.0644124,    1,   3,   0;
      0.00529376,   1,   4,   0;
      0.00891203,   1,   5,   0] "<coeff, power of T/T_red, power of d/d_red, 0>",
    R0=1.02,
    z=0.065,
    xi_0=0.14e-9,
    Gamma_0=0.053,
    qd_inverse=0.4e-9,
    T_ref=608.34) "Coefficients for the thermal conductivity";

  final constant HelmholtzMedia.Interfaces.Types.DynamicViscosityCoefficients
  dynamicViscosityCoefficientsAmmonia(
    dynamicViscosityModel=DynamicViscosityModel.VS1,
    collisionIntegralModel=CollisionIntegralModel.CI1,
    sigma=0.2957 "Lennard-Jones size parameter [nm]",
    epsilon_kappa=386.0 "Lennard-Jones energy parameter [K]",
    a=[
      0.39175,     0;
     -0.59918,     1;
     -0.00022,     2;
      0.19871,     3;
     -0.06942,     4] "collision integral <coeff, power of log(T/(eps/k))>",
    CET=[
      0.0881362,   0.5] "Chapman-Enskog term, 0.021357*sqrt(MM in g/mol)",
    b=[
     -19.572881,    0.0;
      219.73999,   -0.25;
     -1015.3226,   -0.5;
      2471.0125,   -0.75;
     -3375.1717,   -1.0;
      2491.6597,   -1.25;
     -787.26086,   -1.5;
      14.085455,   -2.5;
     -0.34664158,  -5.5] "second viscosity virial coefficient <coeff, power of T/(eps/k)>",
    reducingTemperature_residual=405.56,
    reducingMolarVolume_residual=1/13696 "1/(13.696 mol/l)",
    reducingViscosity_residual=1.0,
    e=[
      0.0393308,    0.5,   2/3,    0,  0;
      16.724735,    0.5,   5/3,    0,  0;
      1.1975934,    0.5,  14/3,    0,  0;
      0.0016995,   -3.5,  26/3,    0,  0;
     -4.2399794,    1.5,   8/3,    0,  0] "residual <coeff, power of T/T_red, power of d/d_red, 0, 0>")
  "Coefficients for the dynamic viscosity";

  final constant HelmholtzMedia.Interfaces.Types.SurfaceTensionCoefficients
  surfaceTensionCoefficientsAmmonia(
    coeffs=[
      0.1028,     1.211;
     -0.09453,    5.585]) "Coefficients for the surface tension <sigma0 and n>";

  final constant HelmholtzMedia.Interfaces.Types.AncillaryCoefficients
  ancillaryCoefficientsAmmonia(
    T_crit=405.56,
    p_crit=11.3634e6,
    d_crit=13.696*17.03052,
    pressureSaturationModel=PressureSaturationModel.PS5,
    pressureSaturation=[
      -7.3128,     1.0;
       3.8888,     1.5;
      -2.9908,     1.6;
      -2.8636,     3.7] "Gao et al. (2023), Eq. (20)",
    densityLiquidModel=DensityLiquidModel.DL1,
    densityLiquid=[
       0.051236,   0.07;
       3.7925,     0.46;
      -3.5929,     0.77;
       4.6409,     1.05;
      -1.9893,     1.25;
       1.5978,     8.0] "Gao et al. (2023), Eq. (21)",
    densityVaporModel=DensityVaporModel.DV3,
    densityVapor=[
      -0.089966,   0.112;
      -3.8722,     0.473;
      -8.1183,     1.5;
     -25.293,      3.875;
     -54.279,      8.0;
    -400.83,      20.0] "Gao et al. (2023), Eq. (22)",
    pressureMeltingModel=PressureMeltingModel.ML1,
    T_reducing=195.49,
    p_reducing=1e6,
    pressureMelting1=[
       0.00605339-516.613763,   0.0;
       516.613763,              4.4] "Gao et al. (2023), Eq. (25), Simon-Glatzel, solid phase I (up to 217.35 K and 307 MPa)",
    pressureMelting2=fill(0.0, 0, 2),
    pressureMelting3=fill(0.0, 0, 2))
  "Coefficients for the ancillary equations (PS5, DL1, DV3, ML1)";

  annotation (Documentation(info="<html>
<p>
These are the coefficients for ammonia (R-717).
The reference state is set to IIR (saturated liquid at 0 &deg;C: h = 200 kJ/kg, s = 1 kJ/(kg K)).
Note that the REFPROP default for ammonia is different (u = 0 and s = 0 for the saturated liquid at the triple point).
</p>
<p>
The equation of state of Gao et al. (2023) contains two associating terms
n*delta^d*tau^t*exp(eta*(delta-epsilon)^2 + 1/(beta*(tau-gamma)^2 + b)), see <code>residualAssociating</code>.
The ancillary equations (vapour pressure, saturated densities, Eqs. (20)-(22)) and the melting curve of solid phase I
(Simon-Glatzel form, Eq. (25), valid up to 217.35 K and 307 MPa) are those of Gao et al. (2023).
</p>

<dl>
<dt> Gao, K., Wu, J., Bell, I.H., Harvey, A.H. and Lemmon, E.W.</dt>
<dd> <b>A Reference Equation of State with an Associating Term for the Thermodynamic Properties of Ammonia</b><br />
     J. Phys. Chem. Ref. Data, 52:013102, 2023.<br />
     DOI: <a href=\"https://doi.org/10.1063/5.0128269\">10.1063/5.0128269</a>
</dd>
<dt> Monogenidou, S.A., Assael, M.J. and Huber, M.L.</dt>
<dd> <b>Reference Correlation for the Viscosity of Ammonia from the Triple Point to 725 K and up to 50 MPa</b><br />
     J. Phys. Chem. Ref. Data, 47:023102, 2018.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.5036724\">10.1063/1.5036724</a>
</dd>
<dt> Monogenidou, S.A., Assael, M.J. and Huber, M.L.</dt>
<dd> <b>Reference Correlation for the Thermal Conductivity of Ammonia from the Triple-Point Temperature to 680 K and Pressures up to 80 MPa</b><br />
     J. Phys. Chem. Ref. Data, 47:043101, 2018.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.5053087\">10.1063/1.5053087</a>
</dd>
<dt> Mulero, A., Cachadi&ntilde;a, I. and Parra, M.I.</dt>
<dd> <b>Recommended Correlations for the Surface Tension of Common Fluids</b><br />
     J. Phys. Chem. Ref. Data, 41:043105, 2012.<br />
     DOI: <a href=\"https://doi.org/10.1063/1.4768782\">10.1063/1.4768782</a>
</dd>
<dt> Lemmon, E.W., Bell, I.H., Huber, M.L. and McLinden, M.O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP, Version 10.0</b><br />
     National Institute of Standards and Technology, 2018.<br />
     DOI: <a href=\"https://doi.org/10.18434/T4/1502528\">10.18434/T4/1502528</a>
</dd>
</dl>
</html>"));
end Ammonia;
