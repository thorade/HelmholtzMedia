within HelmholtzMedia.Examples;
model PumpingSystem "Modelica.Fluid pumping system example with liquid ethanol"
  extends Modelica.Fluid.Examples.PumpingSystem(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol,
    reservoir(height=5),
    PT1(T=5));
  annotation (experiment(StopTime=2000, Interval=0.4, Tolerance=1e-6), Documentation(info="<html>
<p>
The pumping system example of the Modelica Standard Library (Modelica.Fluid.Examples.PumpingSystem) with ethanol
instead of water: a pump lifts the liquid from a source through a 100 m long pipe (50 m height difference) into a
reservoir, controlled by an on/off controller on the reservoir pressure; a user valve opens at t = 200 s.
Ethanol is liquid at ambient pressure and 20 &deg;C, so no fluid parameters need to be changed.
Compared to the original, the reservoir is higher (5 m instead of 3 m) and the first-order filter of the
measured reservoir pressure is slower (T = 5 s); with the original values, the pump switched
on and off too quickly after start-up and the simulation failed at the closed check valve of the pump.
</p>
</html>"));
end PumpingSystem;
