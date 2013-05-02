within HelmholtzMedia.Examples.Validation;
model ReferenceState
  package Medium = HelmholtzMedia.HelmholtzFluids.HMDS;
  Medium.SpecificEnthalpy h_ref;
  Medium.SpecificEntropy s_ref;

  // remember to delete file manually before starting!
protected
  String fileName = "ReferenceState.csv";
  Medium.SaturationProperties sat;
  final constant Medium.Temperature T_IIR = 273.15; // 0°C;
  final constant Medium.Temperature T_ASHRAE = 233.15; // -40°C;
  final constant Medium.AbsolutePressure p_NBP = 101325; // 1.01325 bar = 1 atm

algorithm
  // sat := Medium.setSat_T(T=T_IIR);
  // sat := Medium.setSat_T(T=T_ASHRAE);
   sat := Medium.setSat_p(p=p_NBP);

  s_ref := sat.liq.s;
  h_ref := sat.liq.h;

  // While csv originally stood for comma-seperated-values, MS Excel uses semicolons to seperate the values
  // Modelica.Utilities.Streams.print("idealPower[1,1];" + "s_ref;" + "idealPower[2,1];" + "h_ref;", fileName);
  // Modelica.Utilities.Streams.print("idealPower1;" + "sref;" + "idealPower2;" + "href;", fileName);
  Modelica.Utilities.Streams.print( String(Medium.helmholtzCoefficients.idealPower[1,1],significantDigits=30)+";"
                                   +String(s_ref,significantDigits=30)+";"
                                   +String(Medium.helmholtzCoefficients.idealPower[2,1],significantDigits=30)+";"
                                   +String(h_ref,significantDigits=30)+";",
                                    fileName);

annotation (
Documentation(info="<html>
The choice of reference state is more or less arbitrary,
there are at least three standard reference states:

<dl>

<dt><b>International Institute of Refrigeration (IIR)</b></dt>
<dd>at T=0°C and saturated liquid set h=200 kJ/kg and s=1 kJ/kgK</dd>

<dt><b>American Society of Heating, Refrigerating and Air-Conditioning Engineers (ASHRAE)</b></dt>
<dd>at T=-40°C and saturated liquid set h=0 kJ/kg and s=0 kJ/kgK</dd>

<dt><b>normal boiling point (NBP)</b></dt>
<dd>at p=1 atm = 1.01325 bar  and saturated liquid set h=0 kJ/kg and s=0 kJ/kgK</dd>

</dl>
</html>"),
experiment(NumberOfIntervals=1));
end ReferenceState;
