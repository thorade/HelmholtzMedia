within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function cp0 "returns ideal gas heat capacity from EoS"
  input HelmholtzDerivs f;
  output SpecificHeatCapacity cp0;

algorithm
  cp0 := f.R_s*(1 - f.tau*f.tau*f.itt);
end cp0;
