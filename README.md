## HelmholtzFluids
This library is written in Modelica.
The purpose of this library is to calculate fluid properties from an equation of state (EoS), directly within Modelica and not from an external dll.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.
I am planning to submit this library to the Modelica 2012 conference.

### ToDo
* Improve `BaseProperties` (should independentVariables be used here?)
* Add `annotation(Documentation(info=...html...))` using copy+paste from conference paper
* Add non-analytic critical terms for residual Helmholtz energy (needed for water and CO2)
* Add hyperbolic terms for ideal Helmholtz energy (used in short technical EoS and GERG-2008)
* Add choice of reference state
* Add `setState_hs` (needed for turbine calculation, when power is given and p_out is to be determined)
* And of course: more testing!
  
  
### Implemented Fluids
* Butane (n-Butane)
* Ethanol
* Isobutane
* Isopentane
* Propane
* R134a

### Fluids that might be included in future versions
* Butane (short technical EoS)
* HFC-125 (should be easy)
* HFC-1234yf (ECS viscosity model not yet implemented)
* CO2 (non-analytical terms not yet implemented)
* water according to IAPWS-95 (non-analytical terms not yet implemented)
* nitrogen
* methane
* ammonia
* working fluids for Organic-Rankine-Cycles.