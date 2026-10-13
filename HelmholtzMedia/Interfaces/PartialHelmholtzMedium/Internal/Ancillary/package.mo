within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.Internal;
package Ancillary "Helper functions for determining phase regions and start values"

  annotation (Documentation(info="<html>
<p>
Helper functions for determining the phase region and for computing start values of the iterative functions
(setState_*, setSat_*). They are explicit and cheap, but only approximate:
</p>
<ul>
<li>ancillary equations for the saturation pressure and the saturated liquid and vapour densities, and their inverses,</li>
<li>the melting curve,</li>
<li>start values from the van der Waals and Soave equations of state (density_pT_Soave, pressure_dT_Waals, temperature_pd_Waals),</li>
<li>start values of the saturated densities for the Newton iteration in setSat_T (setSat_T_start).</li>
</ul>
<p>
The coefficients of the ancillary equations are given per fluid in ancillaryCoefficients (record HelmholtzMedia.Interfaces.Types.AncillaryCoefficients).
</p>
</html>"));
end Ancillary;
