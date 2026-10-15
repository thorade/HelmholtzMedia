# Fluids in HelmholtzMedia

All fluids are in [`HelmholtzMedia.HelmholtzFluids`](HelmholtzMedia/HelmholtzFluids).
Each fluid package documents its coefficients and references; this page gives an overview.
Numbers in brackets refer to the [list of references](#references) below.

| Fluid | Reference state | Equation of state | Viscosity | Thermal conductivity | Surface tension |
|---|---|---|---|---|---|
| [Ammonia](HelmholtzMedia/HelmholtzFluids/Ammonia/package.mo) | IIR | [33] | [34] | [35] | [15] |
| [Butane](HelmholtzMedia/HelmholtzFluids/Butane/package.mo) | IIR | [1] | [2] | [3] | [15] |
| [Carbondioxide](HelmholtzMedia/HelmholtzFluids/Carbondioxide/package.mo) | IIR | [4] | [5] | [6] | [15] |
| [Carbondioxide_Short](HelmholtzMedia/HelmholtzFluids/Carbondioxide_Short/package.mo) | IIR | [7] | [5] | [6] | [15] |
| [Ethanol](HelmholtzMedia/HelmholtzFluids/Ethanol/package.mo) | NBP | [8] | [9] | [10] | [11] |
| [Helium](HelmholtzMedia/HelmholtzFluids/Helium/package.mo) | NBP | [12] | [13] | [14] | [15] |
| [Hexamethyldisiloxane](HelmholtzMedia/HelmholtzFluids/Hexamethyldisiloxane/package.mo) | NBP | [16] | – | – | [17] |
| [Isobutane](HelmholtzMedia/HelmholtzFluids/Isobutane/package.mo) | IIR | [1] | [18] | [19] | [15] |
| [Isopentane](HelmholtzMedia/HelmholtzFluids/Isopentane/package.mo) | NBP | [20] | n.d. | n.d. | [15] |
| [Pentane](HelmholtzMedia/HelmholtzFluids/Pentane/package.mo) | NBP | [21] | n.d. | n.d. | [15] |
| [Propane](HelmholtzMedia/HelmholtzFluids/Propane/package.mo) | IIR | [22] | [23] | [24] | [15] |
| [Propylene](HelmholtzMedia/HelmholtzFluids/Propylene/package.mo) | IIR | [36] | – | [37] | [15] |
| [R134a](HelmholtzMedia/HelmholtzFluids/R134a/package.mo) | IIR | [25] | [26] | [27] | [15] |
| [R134a_ASHRAE](HelmholtzMedia/HelmholtzFluids/R134a_ASHRAE/package.mo) | ASHRAE | [25] | [26] | [27] | [15] |
| [R134a_NBP](HelmholtzMedia/HelmholtzFluids/R134a_NBP/package.mo) | NBP | [25] | [26] | [27] | [15] |
| [R32](HelmholtzMedia/HelmholtzFluids/R32/package.mo) | IIR | [28] | [29] | [30] | [15] |

– not available (calling the function raises an error)  
n.d. implemented, but the source of the coefficients is not documented

**Reference states** (enthalpy h and entropy s of the saturated liquid):

- IIR: h = 200 kJ/kg and s = 1 kJ/(kg K) at 0 °C
- ASHRAE: h = 0 and s = 0 at −40 °C
- NBP: h = 0 and s = 0 at the normal boiling point (101.325 kPa)

Where the original publication uses a different reference state, the ideal-gas constants are adjusted so that the REFPROP default reference state of the fluid is met exactly.

**Notes**

- All fluids use the molar gas constant `Modelica.Constants.R` = 8.31446261815324 J/(mol K). Equations fitted with a different value (e.g. Helium, 8.314472 J/(mol K)) therefore deviate by about 1e-6, far below their uncertainty.
- The thermal conductivity correlations include the simplified critical enhancement of Olchowy and Sengers [31], except Helium, which has its own critical term [14].
- Coefficients are taken from the publications or from REFPROP [32].

## References

1. D. Bücker and W. Wagner. Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane. *J. Phys. Chem. Ref. Data*, 35(2):929–1019, 2006. doi:[10.1063/1.1901687](https://doi.org/10.1063/1.1901687)
2. E. Vogel, C. Küchenmeister and E. Bich. Viscosity correlation for n-butane in the fluid region. *High Temp. – High Press.*, 31(2):173–186, 1999. doi:[10.1068/htrt154](https://doi.org/10.1068/htrt154)
3. R. A. Perkins, M. L. V. Ramires, C. A. Nieto de Castro and L. Cusco. Measurement and Correlation of the Thermal Conductivity of Butane from 135 K to 600 K at Pressures to 70 MPa. *J. Chem. Eng. Data*, 47(5):1263–1271, 2002. doi:[10.1021/je0101202](https://doi.org/10.1021/je0101202)
4. R. Span and W. Wagner. A New Equation of State for Carbon Dioxide Covering the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa. *J. Phys. Chem. Ref. Data*, 25(6):1509–1596, 1996. doi:[10.1063/1.555991](https://doi.org/10.1063/1.555991)
5. A. Fenghour, W. A. Wakeham and V. Vesovic. The Viscosity of Carbon Dioxide. *J. Phys. Chem. Ref. Data*, 27(1):31–44, 1998. doi:[10.1063/1.556013](https://doi.org/10.1063/1.556013)
6. V. Vesovic, W. A. Wakeham, G. A. Olchowy, J. V. Sengers, J. T. R. Watson and J. Millat. The Transport Properties of Carbon Dioxide. *J. Phys. Chem. Ref. Data*, 19(3):763–808, 1990. doi:[10.1063/1.555875](https://doi.org/10.1063/1.555875)
7. R. Span and W. Wagner. Equations of State for Technical Applications. III. Results for Polar Fluids. *Int. J. Thermophys.*, 24(1):111–162, 2003. doi:[10.1023/A:1022362231796](https://doi.org/10.1023/A:1022362231796)
8. J. A. Schroeder, S. G. Penoncello and J. S. Schroeder. A Fundamental Equation of State for Ethanol. *J. Phys. Chem. Ref. Data*, 43(4):043102, 2014. doi:[10.1063/1.4895394](https://doi.org/10.1063/1.4895394)
9. S. B. Kiselev, J. F. Ely, I. M. Abdulagatov and M. L. Huber. Generalized SAFT-DFT/DMT Model for the Thermodynamic, Interfacial, and Transport Properties of Associating Fluids: Application for n-Alkanols. *Ind. Eng. Chem. Res.*, 44(17):6916–6927, 2005. doi:[10.1021/ie050010e](https://doi.org/10.1021/ie050010e)
10. M. J. Assael, E. A. Sykioti, M. L. Huber and R. A. Perkins. Reference Correlation of the Thermal Conductivity of Ethanol from the Triple Point to 600 K and up to 245 MPa. *J. Phys. Chem. Ref. Data*, 42(2):023102, 2013. doi:[10.1063/1.4797368](https://doi.org/10.1063/1.4797368)
11. A. Mulero, I. Cachadiña and E. L. Sanjuán. Surface Tension of Alcohols. Data Selection and Recommended Correlations. *J. Phys. Chem. Ref. Data*, 44(3):033104, 2015. doi:[10.1063/1.4927858](https://doi.org/10.1063/1.4927858)
12. D. O. Ortiz Vega, K. R. Hall, J. C. Holste, A. H. Harvey and E. W. Lemmon. An Equation of State for the Thermodynamic Properties of Helium. *NIST Internal Report 8474*, 2023. doi:[10.6028/NIST.IR.8474](https://doi.org/10.6028/NIST.IR.8474)
13. V. D. Arp, R. D. McCarty and D. G. Friend. Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa. *NIST Technical Note 1334 (revised)*, 1998. [PDF](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf)
14. B. A. Hands and V. D. Arp. A Correlation of Thermal Conductivity Data for Helium. *Cryogenics*, 21(12):697–703, 1981. doi:[10.1016/0011-2275(81)90211-3](https://doi.org/10.1016/0011-2275(81)90211-3)
15. A. Mulero, I. Cachadiña and M. I. Parra. Recommended Correlations for the Surface Tension of Common Fluids. *J. Phys. Chem. Ref. Data*, 41(4):043105, 2012. doi:[10.1063/1.4768782](https://doi.org/10.1063/1.4768782)
16. P. Colonna, N. R. Nannan, A. Guardone and E. W. Lemmon. Multiparameter Equations of State for Selected Siloxanes. *Fluid Phase Equilib.*, 244(2):193–211, 2006. doi:[10.1016/j.fluid.2006.04.015](https://doi.org/10.1016/j.fluid.2006.04.015)
17. A. Mulero and I. Cachadiña. Recommended Correlations for the Surface Tension of Several Fluids Included in the REFPROP Program. *J. Phys. Chem. Ref. Data*, 43(2):023104, 2014. doi:[10.1063/1.4878755](https://doi.org/10.1063/1.4878755)
18. E. Vogel, C. Küchenmeister and E. Bich. Viscosity Correlation for Isobutane over Wide Ranges of the Fluid Region. *Int. J. Thermophys.*, 21(2):343–356, 2000. doi:[10.1023/A:1006623310780](https://doi.org/10.1023/A:1006623310780)
19. R. A. Perkins. Measurement and Correlation of the Thermal Conductivity of Isobutane from 114 K to 600 K at Pressures to 70 MPa. *J. Chem. Eng. Data*, 47(5):1272–1279, 2002. doi:[10.1021/je010121u](https://doi.org/10.1021/je010121u)
20. E. W. Lemmon and R. Span. Short Fundamental Equations of State for 20 Industrial Fluids. *J. Chem. Eng. Data*, 51(3):785–850, 2006. doi:[10.1021/je050186n](https://doi.org/10.1021/je050186n)
21. R. Span and W. Wagner. Equations of State for Technical Applications. II. Results for Nonpolar Fluids. *Int. J. Thermophys.*, 24(1):41–109, 2003. doi:[10.1023/A:1022310214958](https://doi.org/10.1023/A:1022310214958)
22. E. W. Lemmon, M. O. McLinden and W. Wagner. Thermodynamic Properties of Propane. III. A Reference Equation of State for Temperatures from the Melting Line to 650 K and Pressures up to 1000 MPa. *J. Chem. Eng. Data*, 54(12):3141–3180, 2009. doi:[10.1021/je900217v](https://doi.org/10.1021/je900217v)
23. E. Vogel, C. Küchenmeister, E. Bich and A. Laesecke. Reference Correlation of the Viscosity of Propane. *J. Phys. Chem. Ref. Data*, 27(5):947–970, 1998. doi:[10.1063/1.556025](https://doi.org/10.1063/1.556025)
24. K. N. Marsh, R. A. Perkins and M. L. V. Ramires. Measurement and Correlation of the Thermal Conductivity of Propane from 86 K to 600 K at Pressures to 70 MPa. *J. Chem. Eng. Data*, 47(4):932–940, 2002. doi:[10.1021/je010001m](https://doi.org/10.1021/je010001m)
25. R. Tillner-Roth and H. D. Baehr. An International Standard Formulation for the Thermodynamic Properties of 1,1,1,2-Tetrafluoroethane (HFC-134a) for Temperatures from 170 K to 455 K and Pressures up to 70 MPa. *J. Phys. Chem. Ref. Data*, 23(5):657–729, 1994. doi:[10.1063/1.555958](https://doi.org/10.1063/1.555958)
26. M. L. Huber, A. Laesecke and R. A. Perkins. Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a. *Ind. Eng. Chem. Res.*, 42(13):3163–3178, 2003. doi:[10.1021/ie0300880](https://doi.org/10.1021/ie0300880)
27. R. A. Perkins, A. Laesecke, J. Howley, M. L. V. Ramires, A. N. Gurova and L. Cusco. Experimental Thermal Conductivity Values for the IUPAC Round-Robin Sample of 1,1,1,2-Tetrafluoroethane (R134a). *NISTIR 6605*, 2000. doi:[10.6028/NIST.IR.6605](https://doi.org/10.6028/NIST.IR.6605)
28. R. Tillner-Roth and A. Yokozeki. An International Standard Equation of State for Difluoromethane (R-32) for Temperatures from the Triple Point at 136.34 K to 435 K and Pressures up to 70 MPa. *J. Phys. Chem. Ref. Data*, 26(6):1273–1328, 1997. doi:[10.1063/1.556002](https://doi.org/10.1063/1.556002)
29. D. Velliadou, K. D. Antoniadis, M. J. Assael and M. L. Huber. Reference Correlation for the Viscosity of Difluoromethane (R-32) from the Triple Point to 425 K and up to 70 MPa. *Int. J. Thermophys.*, 43(8):129, 2022. doi:[10.1007/s10765-022-03050-9](https://doi.org/10.1007/s10765-022-03050-9)
30. R. A. Perkins and M. L. Huber. Thermal conductivity correlation for R-32, unpublished (2005), as implemented in REFPROP 10 [32]; functional form of [24].
31. G. A. Olchowy and J. V. Sengers. A Simplified Representation for the Thermal Conductivity of Fluids in the Critical Region. *Int. J. Thermophys.*, 10(2):417–426, 1989. doi:[10.1007/BF01133538](https://doi.org/10.1007/BF01133538)
32. E. W. Lemmon, I. H. Bell, M. L. Huber and M. O. McLinden. NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties – REFPROP, Version 10.0. National Institute of Standards and Technology, 2018. doi:[10.18434/T4/1502528](https://doi.org/10.18434/T4/1502528)
33. K. Gao, J. Wu, I. H. Bell, A. H. Harvey and E. W. Lemmon. A Reference Equation of State with an Associating Term for the Thermodynamic Properties of Ammonia. *J. Phys. Chem. Ref. Data*, 52(1):013102, 2023. doi:[10.1063/5.0128269](https://doi.org/10.1063/5.0128269)
34. S. A. Monogenidou, M. J. Assael and M. L. Huber. Reference Correlation for the Viscosity of Ammonia from the Triple Point to 725 K and up to 50 MPa. *J. Phys. Chem. Ref. Data*, 47(2):023102, 2018. doi:[10.1063/1.5036724](https://doi.org/10.1063/1.5036724)
35. S. A. Monogenidou, M. J. Assael and M. L. Huber. Reference Correlation for the Thermal Conductivity of Ammonia from the Triple-Point Temperature to 680 K and Pressures up to 80 MPa. *J. Phys. Chem. Ref. Data*, 47(4):043101, 2018. doi:[10.1063/1.5053087](https://doi.org/10.1063/1.5053087)
36. E. W. Lemmon, U. Overhoff, M. O. McLinden and W. Wagner. A Reference Equation of State for the Thermodynamic Properties of Propene for Temperatures from the Melting Line to 575 K and Pressures up to 1000 MPa, to be submitted to *J. Phys. Chem. Ref. Data*; coefficients as implemented in REFPROP 10 [32].
37. M. J. Assael, A. Koutian, M. L. Huber and R. A. Perkins. Reference Correlations of the Thermal Conductivity of Ethene and Propene. *J. Phys. Chem. Ref. Data*, 45(3):033104, 2016. doi:[10.1063/1.4958984](https://doi.org/10.1063/1.4958984). Without viscosity, the critical enhancement is not computed.
