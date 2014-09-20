within HelmholtzMedia.Examples.Validation;
model ReferenceState
  package Medium = HelmholtzMedia.HelmholtzFluids.R134a;

  input String fileName = "ReferenceState_.csv";
  input String separator = ";";
  Medium.ReferenceState ref=Medium.ReferenceState.IIR;
  output Medium.SpecificEnthalpy h_ref;
  output Medium.SpecificEntropy s_ref;

protected
  Boolean hasBeenExecuted;
  Medium.SaturationProperties sat;
  final constant Medium.Temperature T_IIR = 273.15; // 0°C;
  final constant Medium.Temperature T_ASHRAE = 233.15; // -40°C;
  final constant Medium.AbsolutePressure p_NBP = 101325; // 1.01325 bar = 1 atm

algorithm
  if ref == Medium.ReferenceState.IIR then
    sat := Medium.setSat_T(T=T_IIR);
  elseif ref == Medium.ReferenceState.ASHRAE then
    sat := Medium.setSat_T(T=T_ASHRAE);
  elseif ref == Medium.ReferenceState.NBP then
    sat := Medium.setSat_p(p=p_NBP);
  else
    assert(false, "No Ref type selected");
  end if;

  s_ref := Medium.bubbleEntropy(sat);
  h_ref := Medium.bubbleEnthalpy(sat);

  if not hasBeenExecuted then
  if not Modelica.Utilities.Files.exist(fileName) then
    Modelica.Utilities.Streams.print("idealPower1" + separator
                                    +"sref" + separator
                                    +"idealPower2" + separator
                                    +"href" + separator,
                                     fileName);
  end if;
  Modelica.Utilities.Streams.print( String(Medium.helmholtzCoefficients.idealPower[1,1],significantDigits=15) + separator
                                   +String(s_ref,significantDigits=15) + separator
                                   +String(Medium.helmholtzCoefficients.idealPower[2,1],significantDigits=15) + separator
                                   +String(h_ref,significantDigits=15) + separator,
                                    fileName);
  hasBeenExecuted := true;
  end if;

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
experiment(NumberOfIntervals=1, Tolerance=1e-012),
    __Dymola_experimentSetupOutput(doublePrecision=true));
end ReferenceState;
