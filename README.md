# HelmholtzMedia
Modelica library for the calculation of fluid properties from an equation of state (EoS).

## Library description
This library calculates fluid properties from an equation of state (EoS) directly within Modelica. 
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.  

A general description of the library can be found in the related publications:
* Thorade, M. and Saadat, A. (2012). "[HelmholtzMedia - A fluid properties library][1]", 
_Proceedings of the 9th International Modelica Conference_, 
doi:10.3384/ecp1207663
* Thorade, M. (2012) "[Poster HelmholtzMedia][2]", 
_VDI Thermodynamik-Kolloquium 2012_  
* Thorade, M. and Saadat, A. (2013). "[Partial derivatives of thermodynamic state properties for dynamic simulation][3]", 
_Environmental Earth Sciences_, 
doi:10.1007/s12665-013-2394-z

The following fluids have been implemented with EoS and transport properties:  
[Butane](HelmholtzFluids/Butane/package.mo), 
[Ethanol](HelmholtzFluids/Ethanol/package.mo), 
[Isobutane](HelmholtzFluids/Isobutane/package.mo), 
[Isopentane](HelmholtzFluids/Isopentane/package.mo), 
[Pentane](HelmholtzFluids/Pentane/package.mo), 
[Propane](HelmholtzFluids/Propane/package.mo), 
[R134a](HelmholtzFluids/R134a/package.mo) (with three reference states). 

The following fluids have been implemented with EoS only:  
[Helium](HelmholtzFluids/Helium/package.mo) (has convergence problems in the far super-critical region),  [Hexamethyldisiloxane (HMDS)](HelmholtzFluids/HMDS/package.mo) (has convergence problems close to critical point). 

## Current release
Download the newest [tagged version](../../tags).  
In the future, there might be a release branch and official releases.

## License
Copyright © 2009-2013 Helmholtz Centre Potsdam, GFZ German Research Centre for Geosciences

This Modelica package is free software and the use is completely at your own risk;  
it can be redistributed and/or modified under the terms of the [Modelica License 2](https://www.modelica.org/licenses/ModelicaLicense2) or newer.  
For license conditions (including the disclaimer of warranty) visit [https://www.modelica.org/licenses/](https://www.modelica.org/licenses/).

## Development and contribution
Any feedback regarding the library is highly appreciated (beta-testers welcome!).  
You may report feedback, issues or feature-requests using the [Issues](../../issues) button.  
Contributions in the form of [Pull Requests](../../pulls) are always welcome.


[1]: http://goo.gl/Ynuky "Conference Paper: HelmholtzMedia implementation"
[2]: http://goo.gl/HeUzM "Conference Poster: HelmholtzMedia implementation"
[3]: http://goo.gl/HsDXN "ISI Journal Paper: Partial derivatives"