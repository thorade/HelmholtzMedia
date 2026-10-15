# HelmholtzMedia
Modelica library for the calculation of fluid properties from an equation of state (EoS).

## Library description
This library calculates fluid properties from an equation of state (EoS) directly within Modelica.
It supports EoS of the form f=f(T,d) meaning Helmholtz energy as a funtion of temperature and density.
In addition to all state properties, this library calculates viscosity, thermal conductivity and surface tension.  
Unless stated otherwise in the fluid documentation, the reference state for enthalpy and entropy is the REFPROP default of the respective fluid.  

A general description of the library can be found in this [poster][4] or in the related publications:
* Thorade, M. and Saadat, A. (2012). "HelmholtzMedia - A fluid properties library",
_Proceedings of the 9th International Modelica Conference_,
doi:[10.3384/ecp1207663][1]
* Thorade, M. and Saadat, A. (2013). "Partial derivatives of thermodynamic state properties for dynamic simulation",
_Environmental Earth Sciences_,
doi:[10.1007/s12665-013-2394-z][2]
* Thorade, M. (2014). "Entropiebasierte Bewertungskriterien für den Wärmeübergang in Kraftwerksprozessen und ihre Relevanz für praktische Anwendungen", 
_Dissertation (TU Hamburg-Harburg)_,
doi:[10.15480/882.1207][3]

The available fluids, their reference states and the cited references for the equation of state, viscosity, thermal conductivity and surface tension are listed in [Fluids.md](Fluids.md).

## Current release
Download the newest [tagged version](https://github.com/thorade/HelmholtzMedia/tags).  
In the future, there might be a release branch and official releases.

## License
Copyright  2009-2026 Matthis Thorade

This Modelica package is free software and the use is completely at your own risk; 
it is available under the [BSD 3-Clause license](LICENSE.md). 
Upon request, it is also available under other licenses, including the [Modelica License 2](https://www.modelica.org/licenses/ModelicaLicense2).  

## Development and contribution
You may report feedback, issues or feature-requests using the [Issues](https://github.com/thorade/HelmholtzMedia/issues) button.  
Code contributions are very welcome, especially in the form of [Pull Requests](https://github.com/thorade/HelmholtzMedia/pulls).


[1]: http://dx.doi.org/10.3384/ecp1207663 "Modelica Conference Paper: HelmholtzMedia"
[2]: http://dx.doi.org/10.1007/s12665-013-2394-z "ISI Journal Paper: Partial derivatives"
[3]: http://dx.doi.org/10.15480/882.1207 "Dissertation"
[4]: http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:247376:1/component/escidoc:247375/20907.pdf "Poster"
