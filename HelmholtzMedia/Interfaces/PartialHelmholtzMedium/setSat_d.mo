within HelmholtzMedia.Interfaces.PartialHelmholtzMedium;
function setSat_d
  "iterative calculation of saturation properties from EoS with Newton-Raphson algorithm"
  input Density d;
  output SaturationProperties sat;

protected
  constant MolarMass MM = fluidConstants[1].molarMass;
  constant SpecificHeatCapacity R_s=Modelica.Constants.R/MM
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

  Real RES[2] "residual function vector";
  Real Jacobian[2,2] "Jacobian matrix";
  Real NS[2] "Newton step vector";

  constant Real lambda(min=0.1,max=1) = 1 "convergence speed, default=1";
  constant Real tolerance=1e-6 "tolerance for p and T and d, needs to be smaller than simulation tolerance";
  Integer iter = 0;
  constant Integer iter_max = 200;

algorithm
  // Modelica.Utilities.Streams.print(" ", "printlog.txt");
  // Modelica.Utilities.Streams.print("setSat_d: d="+String(d),"printlog.txt");

  sat.Tsat := Ancillary.saturationTemperature_d(d=d);
  if (d<d_crit) and (d>dv_trip) then
    // Modelica.Utilities.Streams.print("d<d_crit: input is on vapour side: find d_liq and T_sat", "printlog.txt");
    sat.liq.d := 1.02*Ancillary.bubbleDensity_T(sat.Tsat);
    sat.vap.d := d; // d_vap is a constant

    // calculate residuals: liq-vap (=var-const)
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    RES := {EoS.p(fl)-EoS.p(fv), EoS.g(fl)-EoS.g(fv)};

  while (iter<iter_max) and (iter<1 or abs(RES[1])>tolerance or abs(RES[2])>tolerance or abs(NS[1])>tolerance or abs(NS[2])>tolerance) loop
      iter := iter+1;

      // calculate Jacobian matrix and Newton Step vector
      Jacobian := [EoS.dpdT(fl), EoS.dpTd(fl)-EoS.dpTd(fv);
                   EoS.dgdT(fl), EoS.dgTd(fl)-EoS.dgTd(fv)];
      NS := -Modelica.Math.Matrices.solve(Jacobian,RES);

      // calculate better values for sat.liq.d and sat.Tsat
      sat.liq.d := sat.liq.d +lambda*NS[1];
      sat.Tsat  := sat.Tsat  +lambda*NS[2];

      // check bounds
      sat.liq.d := max(sat.liq.d, 0.98*d_crit);
      sat.liq.d := min(sat.liq.d, 1.02*dl_trip);
      sat.Tsat  := max(sat.Tsat,  0.98*T_trip);
      sat.Tsat  := min(sat.Tsat,  1.02*T_crit);

      // calculate new residuals: liq-vap
      fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
      fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
      RES := {EoS.p(fl)-EoS.p(fv), EoS.g(fl)-EoS.g(fv)};
    end while;
    sat.liq  := setState_dTX(d=sat.liq.d, T=sat.Tsat, phase=1);
    sat.vap  := setState_dTX(d=sat.vap.d, T=sat.Tsat, phase=1);
    sat.psat := sat.liq.p;

  elseif (d>d_crit) and (d<dl_trip) then
    // Modelica.Utilities.Streams.print("d>d_crit: input is on liquid side: find d_vap and T_sat", "printlog.txt");
    sat.vap.d := 0.98*Ancillary.dewDensity_T(sat.Tsat);
    sat.liq.d := d; // d_liq is a constant

    // calculate residuals: vap-liq (=var-const)
    fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
    fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
    RES := {EoS.p(fv)-EoS.p(fl), EoS.g(fv)-EoS.g(fl)};

  while (iter<iter_max) and (iter<1 or abs(RES[1])>tolerance or abs(RES[2])>tolerance or abs(NS[1])>tolerance or abs(NS[2])>tolerance) loop
      iter := iter+1;

      // calculate Jacobian matrix and Newton Step vector
      Jacobian := [EoS.dpdT(fv), EoS.dpTd(fv)-EoS.dpTd(fl);
                   EoS.dgdT(fv), EoS.dgTd(fv)-EoS.dgTd(fl)];
      NS := -Modelica.Math.Matrices.solve(Jacobian,RES);

      // calculate better values for sat.vap.d and sat.Tsat
      sat.vap.d := sat.vap.d +lambda*NS[1];
      sat.Tsat  := sat.Tsat  +lambda*NS[2];

      // check bounds
      sat.vap.d := max(sat.vap.d, 0.98*dv_trip);
      sat.vap.d := min(sat.vap.d, 1.02*d_crit);
      sat.Tsat  := max(sat.Tsat,  0.98*T_trip);
      sat.Tsat  := min(sat.Tsat,  1.02*T_crit);

      // calculate new residuals: vap-liq
      fv := EoS.setHelmholtzDerivsSecond(d=sat.vap.d, T=sat.Tsat, phase=1);
      fl := EoS.setHelmholtzDerivsSecond(d=sat.liq.d, T=sat.Tsat, phase=1);
      RES := {EoS.p(fv)-EoS.p(fl), EoS.g(fv)-EoS.g(fl)};
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
  assert(iter<iter_max, "setSat_d did not converge, input was d=" + String(d)+
                        "; the remaining residuals are RES[1]=" + String(RES[1]) +
                        " and RES[2]=" + String(RES[2]));

end setSat_d;
