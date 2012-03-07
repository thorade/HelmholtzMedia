within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
record SurfaceTensionCoefficients
  Real[:,2] coeffs = fill(0.0, 0, 2) "sigma0 and n";
end SurfaceTensionCoefficients;
