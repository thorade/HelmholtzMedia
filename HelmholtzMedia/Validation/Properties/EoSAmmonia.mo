within HelmholtzMedia.Validation.Properties;
model EoSAmmonia "Ammonia EoS (Gao et al. 2023) against check values computed with CoolProp 8.0.0 (IIR reference state)"
  extends EoSCheck(
    redeclare package Medium = HelmholtzMedia.HelmholtzFluids.Ammonia,
    T={250, 300, 300, 400, 406, 500, 700},
    d={690, 5, 620, 30, 233.25, 100, 1},
    p={42959563.89, 678484.5051, 24208451.62, 4770109.736, 11452452.18, 17436336.08, 341315.079},
    cv={2914.414626, 1927.165827, 2786.407284, 2220.971695, 6050.946721, 2429.233598, 2338.97084},
    w={1846.561714, 417.4763959, 1474.492171, 445.3983377, 225.4457638, 482.1758384, 642.4582718},
    h={130135.0844, 1513242.845, 337881.7739, 1639027.696, 1105261.458, 1723506.941, 2537586.026},
    s={488.7977259, 5591.412333, 1349.583661, 5120.944055, 3518.460495, 4828.880677, 8072.972414},
    Tsat={200, 240, 273.15, 300, 350, 400},
    psat={8609.776626, 102171.0302, 429247.7987, 1061121.502, 3865199.174, 10297199.4},
    dl={728.6693135, 681.4309219, 638.6382049, 600.1699235, 512.4162065, 344.0076036},
    dv={0.08867285906, 0.8969191288, 3.456010881, 8.24427316, 31.33643153, 130.8857835},
    hl={-126220.2624, 49471.08258, 200000, 326380.8388, 582657.1726, 942918.6633},
    hv={1351536.745, 1418647.623, 1461773.198, 1484432.155, 1478168.206, 1289109.072},
    sv={7003.620529, 6119.525895, 5619.341747, 5297.807952, 4769.385564, 3991.508427});
  annotation (experiment(StopTime=1), Documentation(info="<html>
<p>The check values were computed with CoolProp 8.0.0 (equation of state of Gao et al., including the associating terms),
with the IIR reference state. They test in particular the associating terms of the residual Helmholtz energy.</p>
</html>"));
end EoSAmmonia;
