within HelmholtzMedia.Interfaces.PartialHelmholtzMedium.EoS;
function dgTd "returns Gibbs energy derivative (dg/dT)@d=const"
  input EoS.HelmholtzDerivs f;
  output DerEnergyByTemperature dgTd;

algorithm
  dgTd := f.R*(-f.tau*(f.it+f.rt) +(f.i+f.r) +(1+f.delta*f.rd-f.delta*f.tau*f.rtd));
end dgTd;
