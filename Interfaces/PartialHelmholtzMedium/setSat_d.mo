within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setSat_d
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
  input Density d;
  output SaturationProperties sat;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R=Modelica.Constants.R/MM
    "specific gas constant";
  constant Density d_crit=MM/fluidConstants[1].criticalMolarVolume;
  constant Temperature T_trip=fluidConstants[1].triplePointTemperature;
  constant Temperature T_crit=fluidConstants[1].criticalTemperature;
  constant AbsolutePressure p_trip=fluidConstants[1].triplePointPressure;
  constant AbsolutePressure p_crit=fluidConstants[1].criticalPressure;
  constant Density dv_trip = Ancillary.dewDensity_T(T_trip);
  constant Density dl_trip = Ancillary.bubbleDensity_T(T_trip);

  EoS.HelmholtzDerivs fl;
  EoS.HelmholtzDerivs fv;

  AbsolutePressure RES_p;
  SpecificEnergy RES_g;
  DerPressureByDensity dpdT "(dp/dd)@T=const";
  DerPressureByTemperature dpTd "(dp/dT)@d=const";
  DerEnergyByDensity dgdT "(dg/dd)@T=const";
  DerEnergyByTemperature dgTd "(dg/dT)@d=const";
  Real det "determinant of Jacobi matrix";
  constant Real gamma(min=0,max=1) = 1 "convergence speed, default=1";
  constant Real tolerance=1e-6 "tolerance for relative RES_p or RES_g";
  Integer iter = 0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print(" ", "printlog.txt");
  // Modelica.Utilities.Streams.print("setSat_d: d="+String(d),"printlog.txt");

  sat.Tsat := 0.99*Ancillary.saturationTemperature_d(d=d);
  if (d<d_crit) and (d>dv_trip) then
    // Modelica.Utilities.Streams.print("d<d_crit: input is on vapour side: find d_liq and T_sat", "printlog.txt");
    sat.liq.d := Ancillary.bubbleDensity_T(sat.Tsat);
    sat.vap.d := d; // d_vap is a constant

    // calculate residuals: liq-vap (=var-const)
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    RES_p  := EoS.p(fl) - EoS.p(fv);
    RES_g  := EoS.g(fl) - EoS.g(fv);

    while (abs(RES_p/(fv.d*fv.T*fv.R))>tolerance or abs(RES_g/(fv.T*fv.R))>tolerance) and (iter<iter_max) loop
      iter := iter+1;
      // gamma := (iter_max-iter)/iter_max;

      // calculate gradients of residual functions regarding d_liq and T
      dpdT := EoS.dpdT(fl);
      dpTd := EoS.dpTd(fl) - EoS.dpTd(fv);
      dgdT := EoS.dgdT(fl);
      dgTd := EoS.dgTd(fl) - EoS.dgTd(fv);

      // calculate determinant of Jacobi matrix det=ad-bc
      det := dpdT*dgTd-dpTd*dgdT;

      // print for Newton debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("sat.liq.d=" + String(sat.liq.d) + " and sat.Tsat=" + String(sat.Tsat), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and dpdT=" + String(dpdT) + " and dpTd=" + String(dpTd), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_g=" + String(RES_g) + " and dgdT=" + String(dgdT) + " and dgTd=" + String(dgTd), "printlog.txt");
      // Modelica.Utilities.Streams.print("Jacobi determinant det=" +String(det), "printlog.txt");

      // calculate better values for sat.liq.d and sat.Tsat
      sat.liq.d := sat.liq.d -gamma/det*(+dgTd*RES_p -dpTd*RES_g);
      sat.Tsat  := sat.Tsat  -gamma/det*(-dgdT*RES_p +dpdT*RES_g);

      // check bounds
      sat.liq.d := max(sat.liq.d, 0.98*d_crit);
      sat.liq.d := min(sat.liq.d, 1.02*dl_trip);
      sat.Tsat  := max(sat.Tsat,  0.98*T_trip);
      sat.Tsat  := min(sat.Tsat,  1.02*T_crit);

      // calculate new residuals: liq-vap
      fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
      RES_p  := EoS.p(fl) - EoS.p(fv);
      RES_g  := EoS.g(fl) - EoS.g(fv);
    end while;
    sat.liq  := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap  := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);
    sat.psat := sat.liq.p;

  elseif (d>d_crit) and (d<dl_trip) then
    // Modelica.Utilities.Streams.print("d>d_crit: input is on liquid side: find d_vap and T_sat", "printlog.txt");
    sat.vap.d := Ancillary.dewDensity_T(sat.Tsat);
    sat.liq.d := d; // d_liq is a constant

    // calculate residuals: vap-liq (=var-const)
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    RES_p  := EoS.p(fv) - EoS.p(fl);
    RES_g  := EoS.g(fv) - EoS.g(fl);

    while (abs(RES_p/(fl.d*fl.T*fl.R))>tolerance or abs(RES_g/(fl.T*fl.R))>tolerance) and (iter<iter_max) loop
      iter := iter+1;

      // calculate gradients of residual functions regarding d_vap and T
      dpdT := EoS.dpdT(fv);
      dpTd := EoS.dpTd(fv) - EoS.dpTd(fl);
      dgdT := EoS.dgdT(fv);
      dgTd := EoS.dgTd(fv) - EoS.dgTd(fl);

      // calculate determinant of Jacobi matrix det=ad-bc
      det := dpdT*dgTd-dpTd*dgdT;

      // print for Newton debugging
      // Modelica.Utilities.Streams.print("Iteration step " +String(iter), "printlog.txt");
      // Modelica.Utilities.Streams.print("sat.vap.d=" + String(sat.vap.d) + " and sat.Tsat=" + String(sat.Tsat), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_p=" + String(RES_p) + " and dpdT=" + String(dpdT) + " and dpTd=" + String(dpTd), "printlog.txt");
      // Modelica.Utilities.Streams.print("RES_g=" + String(RES_g) + " and dgdT=" + String(dgdT) + " and dgTd=" + String(dgTd), "printlog.txt");
      // Modelica.Utilities.Streams.print("Jacobi determinant det=" +String(det), "printlog.txt");

      // calculate better values for sat.vap.d and sat.Tsat
      sat.vap.d := sat.vap.d -gamma/det*(+dgTd*RES_p -dpTd*RES_g);
      sat.Tsat  := sat.Tsat  -gamma/det*(-dgdT*RES_p +dpdT*RES_g);

      // check bounds
      sat.vap.d := max(sat.vap.d, 0.98*dv_trip);
      sat.vap.d := min(sat.vap.d, 1.02*d_crit);
      sat.Tsat  := max(sat.Tsat,  0.98*T_trip);
      sat.Tsat  := min(sat.Tsat,  1.02*T_crit);

      // calculate new residuals: vap-liq
      fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
      fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
      RES_p  := EoS.p(fv) - EoS.p(fl);
      RES_g  := EoS.g(fv) - EoS.g(fl);
    end while;
    sat.liq  := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap  := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);
    sat.psat := sat.liq.p;

  elseif (d>=dl_trip) or (d<=dv_trip) then
    // Modelica.Utilities.Streams.print("d out of two-phase range, return triple point values", "printlog.txt");
    sat.Tsat:= T_trip;
    sat.psat:= p_trip;
    sat.liq := setState_dTX(d=dl_trip, T=T_trip, phase=1);
    sat.vap := setState_dTX(d=dv_trip, T=T_trip, phase=1);
  else
    // Modelica.Utilities.Streams.print("d=d_crit: return critical values", "printlog.txt");
    sat.Tsat := T_crit;
    sat.psat := p_crit;
    sat.liq := setState_dTX(d=d_crit, T=T_crit, phase=1);
    sat.vap := setState_dTX(d=d_crit, T=T_crit, phase=1);
  end if;
  // Modelica.Utilities.Streams.print("setSat_d total iteration steps " + String(iter), "printlog.txt");
  assert(iter<iter_max, "setSat_d did not converge, input was d=" + String(d));

end setSat_d;
