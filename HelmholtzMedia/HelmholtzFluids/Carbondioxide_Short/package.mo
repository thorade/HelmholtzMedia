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
    hasCriticalData=true,
       criticalTemperature=304.307293794,
       criticalPressure=7400132.63572,
       criticalMolarVolume=0.0440098/446.620950521,
       HCRIT0=339320.335796829,
       SCRIT0=1456.67142727164,
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
    T_reducing=304.1282,
    d_reducing=467.6,
    idealLog=[
      +2.50000,         1.],
    idealPower=[
      -6.123273733681842,         0;
       5.114653226097055,         1],
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
Transport properties and surface tension are the same as in Carbondioxide.
The reference state is set to IIR (saturated liquid at 0 &deg;C: h = 200 kJ/kg, s = 1 kJ/(kg K)), the REFPROP default.

<dl>
<dt> Span, R. and Wagner, W.</dt>
<dd> <b>Equations of State for Technical Applications. III. Results for Polar Fluids</b><br />
International Journal of Thermophysics, 24, 111--162 (2003)<br />
     DOI: <a href=\"http://dx.doi.org/10.1023/A:1022362231796\">10.1023/A:1022362231796</a>
</dd>
</dl>
</html>"));

end Carbondioxide_Short;
