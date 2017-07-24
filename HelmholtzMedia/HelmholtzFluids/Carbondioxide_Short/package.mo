within HelmholtzMedia.HelmholtzFluids;
package Carbondioxide_Short "Carbondioxide short FES"
  extends HelmholtzMedia.HelmholtzFluids.Carbondioxide(
    fluidConstants={fluidConstantsCarbondioxide_Short},
    helmholtzCoefficients=helmholtzCoefficientsCarbondioxide_Short);

  final constant FluidConstants
  fluidConstantsCarbondioxide_Short(
    casRegistryNumber="124-38-9",
    iupacName="Carbondioxide" "full name",
    structureFormula="CO2",
    chemicalFormula="CO2",
    molarMass=0.0440098,
    gasConstant=8.31451,
    hasCriticalData=true,
       criticalTemperature=304.1282,
       criticalPressure=7377300,
       criticalMolarVolume=0.0440098/467.6,
       HCRIT0=333885.585791306,
       SCRIT0=1439.02775748364,
    hasIdealGasHeatCapacity=false,
    hasDipoleMoment=true,
       dipoleMoment=0.0,
    hasFundamentalEquation=true,
    hasLiquidHeatCapacity=true,
    hasSolidHeatCapacity=false,
    hasAccurateViscosityData=true,
    hasAccurateConductivityData=true,
    hasVapourPressureCurve=true,
    hasAcentricFactor=true,
       acentricFactor=0.22394,
    triplePointTemperature=216.592,
    triplePointPressure=0.51795e6,
    normalBoilingPoint=194.686,
    meltingPoint=216.592) "Fluid Constants";

  final constant EoS.HelmholtzCoefficients
  helmholtzCoefficientsCarbondioxide_Short(
    idealLog=[
      +2.50000,         1.],
    idealPower=[
      -6.123243569570599,         0;
       5.114633389652261,         1],
    idealEinstein=[
      +1.994270420,       -3.15163;
      +0.621052475,       -6.11190;
      +0.411952928,       -6.77708;
      +1.040289220,      -11.32384;
      +0.0832767753,     -27.08792],
    residualPoly=[
       0.898751080000E+00,  0.25,    1.0;
      -0.212819850000E+01,  1.25,    1.0;
      -0.681903200000E-01,  1.5,     1.0;
       0.763553060000E-01,  0.25,    3.0;
       0.220532530000E-03,  0.875,   7.0],
    residualBwr=[
      +0.415418230000E+00,  2.375,   1.0,     1;
       0.713356570000E+00,  2.0,     2.0,     1;
       0.303542340000E-03,  2.125,   5.0,     1;
      -0.366431430000E+00,  3.5,     1.0,     2;
      -0.144077810000E-02,  6.5,     1.0,     2;
      -0.891667070000E-01,  4.75,    4.0,     2;
      -0.236998870000E-01, 12.5,     2.0,     3],
     residualGauss=fill(0.0, 0, 9),
     residualNonAnalytical=fill(0.0, 0, 12))
  "Coefficients of the Helmholtz EoS";

  annotation (Documentation(info="<html>
  These are the coefficients for Carbondioxide, 
  using a shorter and faster but less accurate equation of state.

<dl>
<dt> Span, R. and Wagner, W.</dt>
<dd> <b>Equations of State for Technical Applications. III. Results for Polar Fluids</b><br />
International Journal of Thermophysics, 24, 111--162 (2003)<br />
     DOI: <a href=\"http://dx.doi.org/10.1023/A:1022362231796\">10.1023/A:1022362231796</a>
</dd>
</dl>
</html>"));

end Carbondioxide_Short;
