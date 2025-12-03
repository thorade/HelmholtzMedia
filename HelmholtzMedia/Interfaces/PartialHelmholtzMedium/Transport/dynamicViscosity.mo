within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Transport;
function dynamicViscosity "Returns dynamic Viscosity"
  input ThermodynamicState state;
  output DynamicViscosity eta;

protected
  DynamicViscosityModel dynamicViscosityModel=dynamicViscosityCoefficients.dynamicViscosityModel;
  constant Real micro=1e-6;

algorithm
  // assert(state.phase <> 2, "dynamicViscosity warning: property not defined in two-phase region", level=AssertionLevel.warning);

  if (dynamicViscosityModel == DynamicViscosityModel.VS0) then
    // hardcoded models return full viscosity in one equation
    eta := dynamicViscosity_residual(state);
  else
    // composite models
    eta := dynamicViscosity_dilute(state) + dynamicViscosity_initial(state) + dynamicViscosity_residual(state);
  end if;
  // RefProp results are in µPa*s where µ means micro or 1E-6 but SI default is Pa*s
  eta := micro*eta;

  annotation (Documentation(info="<html>
<p>
This model is identical to the RefProp VS1 or VS2 model.

The viscosity is split into three contributions:
zero density (dilute gas) viscosity eta_0,
initial density contribution eta_1
and residual contribution eta_r.

This allows to develop functions for each contribution separately.
The so called background viscosity is the sum of initial and residual viscosity.

At the critical point and a small region around the critical point, the viscosity is enhanced.
As this critical enhancement is small, it is neglected here.

Special thanks go to Eric W. Lemmon for answering all my emails
and programming a special version of RefProp that outputs also intermediate values.

</p>

<dl>
<dt> Lemmon, Eric W.; Huber, M. L. and McLinden, M. O.</dt>
<dd> <b>NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties - REFPROP. 9.0</b><br />
     National Institute of Standards and Technology, Standard Reference Data Program. Gaithersburg<br />
     URL: <a href=\"http://www.nist.gov/srd/nist23.cfm\">http://www.nist.gov/srd/nist23.cfm</a>
</dd>
<dt>Vogel, E.; K&uuml;chenmeister, C. and Birch, E.</dt>
<dd> <b>Reference correlation of the viscosity of propane</b>.<br />
     Journal of Thermophysics (1998) 10, 417-426.<br />
     DOI: <a href=\"http://dx.doi.org/10.1007/BF01133538\">10.1007/BF01133538</a>
</dd>
</dl>
</html>"));
end dynamicViscosity;
