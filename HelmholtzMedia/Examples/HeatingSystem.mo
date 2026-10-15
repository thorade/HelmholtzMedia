within HelmholtzMedia.Examples;
model HeatingSystem "Modelica.Fluid heating system example with liquid ethanol at 3 bar"
  extends Modelica.Fluid.Examples.HeatingSystem(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ethanol,
    system(p_ambient=300000),
    tank(ports(each p(start=3.1e5))),
    pump(
      p_a_start=310000,
      p_b_start=330000,
      p_a_nominal=310000,
      p_b_nominal=330000),
    heater(p_a_start=330000),
    radiator(
      p_a_start=310000,
      state_a(p(start=310000)),
      state_b(p(start=310000))),
    pipe(p_a_start=330000));
  annotation (experiment(StopTime=6000, Interval=2, Tolerance=1e-6), Documentation(info="<html>
<p>
The heating system example of the Modelica Standard Library (Modelica.Fluid.Examples.HeatingSystem) with ethanol
instead of water: a closed loop with an open tank, a pump, a heater, a pipe, a valve and a radiator, with heat
transfer to the ambient. The valve opens further at t = 2000 s.
</p>
<p>
At 1 bar, ethanol boils at about 78 &deg;C, while the heater and the pipe start at 80 &deg;C.
To keep the ethanol liquid, the ambient pressure (and therefore the pressure in the open tank) is raised
to 3 bar, where ethanol boils at about 112 &deg;C. All pressure start values and the nominal pressures of the pump
are shifted by +2 bar accordingly; the pump keeps its original pressure rise of 0.2 bar.
</p>
</html>"));
end HeatingSystem;
