within HelmholtzMedia.Examples.Validation;
model idealGasLimit
  replaceable package Medium = HelmholtzFluids.Helium;

  parameter Medium.Temperature T=298.15;
  Modelica.SIunits.SpecificVolume v;
  Medium.Density d=1/v;

  Medium.ThermodynamicState state = Medium.setState_dT(d=d, T=T, phase=1);
  Medium.EoS.HelmholtzDerivs f = Medium.EoS.setHelmholtzDerivsThird(d=d, T=T, phase=1);
//Medium.EoS.HelmholtzDerivs f0 = Medium.EoS.setHelmholtzDerivsThird(d=0, T=T, phase=1);

  constant Medium.MolarMass MM = Medium.fluidConstants[1].molarMass;
  constant Medium.SpecificHeatCapacity R=fluidConstants[1].gasConstant/MM
    "specific gas constant";
  Real Z = (state.p*v)/(R*state.T);

  Medium.SpecificHeatCapacity cp = Medium.specificHeatCapacityCp(state);
  Medium.SpecificHeatCapacity cp0 = Medium.EoS.cp0(f);

//Real B=f.rd*f.delta/d;
  Real B=f.rd/f.d_crit;
//Real B0=f0.rd/f0.d_crit;

//Real C=f.rdd*f.delta^2/d^2;
  Real C=f.rdd/f.d_crit^2;
//Real C0=f0.rdd/f0.d_crit^2;

//Real D=f.rddd*f.delta^3/d^3;
  Real D=f.rddd/f.d_crit^3;
//Real D0=f0.rddd/f0.d_crit^3;

equation
  v = exp(time/10);

  annotation (Documentation(info="
<html>
<style type=\"text/css\">
  code{background-color:#EEEEEE; padding:2px; margin:2px;}
</style>
This model is used to validate the behavior of the EoS in the ideal gas limit, that is at <code> lim d &rarr; 0 </code>.
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

</html>"),
experiment(StopTime=100000));
end idealGasLimit;
