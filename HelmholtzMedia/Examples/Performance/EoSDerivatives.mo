within HelmholtzMedia.Examples.Performance;
model EoSDerivatives "Performance test: repeated evaluation of the Helmholtz energy derivatives"
  extends Modelica.Icons.Example;
  replaceable package Medium = HelmholtzMedia.HelmholtzFluids.Propane;
  parameter Integer level=2 "order of the derivatives passed to setHelmholtzDerivs (0 to 3)";
  parameter Integer nT=50 "number of temperatures";
  parameter Integer nd=50 "number of densities";
  parameter Integer nRepeat=20 "number of repetitions of the whole grid";
  parameter String fileName="perf_EoSDerivatives.txt";
  Real checksum "sum of all computed derivatives";
  Integer nCalls "number of EoS evaluations";

protected
  constant Medium.Temperature T_crit=Medium.fluidConstants[1].criticalTemperature;
  constant Medium.Density d_crit=Medium.fluidConstants[1].molarMass/Medium.fluidConstants[1].criticalMolarVolume;
  constant Medium.Temperature T_min=max(Medium.fluidLimits.TMIN, 0.5*T_crit);
  constant Medium.Temperature T_max=min(Medium.fluidLimits.TMAX, 3*T_crit);
  constant Medium.Density d_min=1e-3*d_crit;
  constant Medium.Density d_max=min(Medium.fluidLimits.DMAX, 3*d_crit);

  function run "the loops are inside a function, since the HelmholtzDerivs record has binding equations"
    input Integer level;
    input Integer nT;
    input Integer nd;
    input Integer nRepeat;
    output Real checksum=0;
    output Integer nCalls=0;
  protected
    Medium.Temperature T;
    Medium.Density d;
    Medium.EoS.HelmholtzDerivs f;
  algorithm
    for k in 1:nRepeat loop
      for i in 1:nT loop
        T := T_min*(T_max/T_min)^((i - 1)/(nT - 1));
        for j in 1:nd loop
          d := d_min*(d_max/d_min)^((j - 1)/(nd - 1));
          if level == 0 then
            f := Medium.EoS.setHelmholtzDerivs(d=d, T=T, order=0);
            checksum := checksum + f.i + f.r;
          elseif level == 1 then
            f := Medium.EoS.setHelmholtzDerivs(d=d, T=T, order=1);
            checksum := checksum + f.i + f.it + f.r + f.rt + f.rd;
          elseif level == 2 then
            f := Medium.EoS.setHelmholtzDerivs(d=d, T=T, order=2);
            checksum := checksum + f.i + f.it + f.itt + f.r + f.rt + f.rtt + f.rtd + f.rd + f.rdd;
          else
            f := Medium.EoS.setHelmholtzDerivs(d=d, T=T, order=3);
            checksum := checksum + f.i + f.it + f.itt + f.ittt + f.r + f.rt + f.rtt + f.rttt + f.rtd + f.rttd + f.rtdd + f.rd + f.rdd + f.rddd;
          end if;
          nCalls := nCalls + 1;
        end for;
      end for;
    end for;
  end run;

algorithm
  when initial() then
    (checksum, nCalls) := run(level, nT, nd, nRepeat);
    Modelica.Utilities.Files.remove(fileName);
    Modelica.Utilities.Streams.print("level=" + String(level) + " nCalls=" + String(nCalls)
      + " checksum=" + String(checksum, significantDigits=17), fileName);
  end when;

  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>
Evaluates the Helmholtz energy and its derivatives (selected by <code>level</code>) on a grid of
temperatures from max(TMIN, 0.5*T_crit) to min(TMAX, 3*T_crit) and densities from 1e-3*d_crit to min(DMAX, 3*d_crit),
both logarithmically spaced, nRepeat times. The EoS is evaluated regardless of phase (no phase check),
which is fine for a performance test. The checksum allows to verify that a different implementation gives the same results.
</p>
</html>"));
end EoSDerivatives;
