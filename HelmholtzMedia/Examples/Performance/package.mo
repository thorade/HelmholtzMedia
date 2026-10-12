within HelmholtzMedia.Examples;
package Performance "Performance tests: repeated evaluation of EoS and property functions"
  extends Modelica.Icons.ExamplesPackage;

  annotation (Documentation(info="<html>
<p>
These models call EoS or property functions many times inside a <code>when initial()</code> algorithm,
so the computational effort shows up as CPU time for initialization (dslog.txt in Dymola).
Each model also writes a checksum of all results to a file, so that a faster implementation
can be checked to give the same results.
</p>
</html>"));
end Performance;
