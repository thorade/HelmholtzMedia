## HelmholtzMedia
This library is written in Modelica.
The purpose of this library is to calculate fluid properties from an equation of state (EoS), directly within Modelica and not from an external dll.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.  
Also see this [overview over the implementation][1].


### Implemented Fluids
* Butane
* Ethanol
* Isobutane
* Isopentane
* Pentane (hyperbolic terms by Martin Ryhl Kærn)
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
* helium
* working fluids for Organic-Rankine-Cycles.

[1]: http://goo.gl/HeUzM "HelmholtzMedia CheatSheet"