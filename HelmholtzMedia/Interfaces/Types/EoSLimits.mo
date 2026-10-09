within HelmholtzMedia.Interfaces.Types;
record EoSLimits "Validity limits of the Helmholtz EoS"
  extends Modelica.Media.Interfaces.Types.FluidLimits(
    PMAX(max=Modelica.Constants.inf));
  annotation (Documentation(info="<html>
<p>Same as <code>Modelica.Media.Interfaces.Types.FluidLimits</code>, but without the upper limit of 100 MPa for <code>PMAX</code>,
because several equations of state are valid up to 1000 MPa.</p>
</html>"));
end EoSLimits;
