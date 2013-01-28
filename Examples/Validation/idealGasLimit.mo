within HelmholtzMedia.Examples.Validation;
model idealGasLimit
  package Medium = HelmholtzFluids.Propane;

  parameter Medium.Temperature T=298.15;
  Medium.Density d;

  Medium.EoS.HelmholtzDerivs f = Medium.EoS.setHelmholtzDerivsThird(d=d, T=T, phase=1);
  Medium.ThermodynamicState state = Medium.setState_dT(d=d, T=T, phase=1);

  Real B=f.rd*f.delta/d;
  Real C=f.rdd*f.delta^2/d^2;
  Real D=f.rddd*f.delta^3/d^3;

equation
  d = 100/exp(time);

  annotation (Documentation(info="
<html>
<style type=\"text/css\">
  code{background-color:#EEEEEE; padding:2px; margin:2px;}
</style>
<body>
This model is used to validate the behavior of the EoS in the ideal gas limit, that is at <code> lim d -> 0 </code>.
<br />
How to use:
<ol>
  <li>Simulate with any <code>T_trip &lt; T &lt; T_max</code></li>
  <li>Plot <code>f.rd</code> and <code>B</code> </li>
  <li>Make <code>d</code> the independent variable of the plot </li>
  <li>Set the horizontal axis <code>(d)</code> to logarithmic scale</li>
  <li>Set line style to dotted and marker style to cross</li>
  <li>Look at the graphs </li>
</ol>

</body>
</html>"),experiment(StopTime=1000, __Dymola_NumberOfIntervals=10000),
    __Dymola_experimentSetupOutput);
end idealGasLimit;
