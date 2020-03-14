## Description
WInc3D provides an integrated wind farm simulation framework that allows detailed analyses of wake–to–wake and turbine–to–wake interactions. The code is built upon the high–order finite difference
numerical solver incompact3D (Laizet and Lamballais, 2009) which makes use of an efficient 2D domain decomposition algorithm (Laizet and Li, 2011) that allows the code to scale on up to O(10$^5$)
computational cores. WInc3D offers a number of built-in models, including a native actuator line model (ALM). The models have been validated using both wind tunnel data (Deskos et al 2019) and supervisory control and data acquisition (SCADA) measurements from a utility–scale offshore wind farm (Deskos et al 2018)
### This repository 


To compile
-----------
In order to compile you need to have mpif90 installed in your computer. 
Once you have it installed just do:
```bash
make 
```

How to cite WInc3D ?
--------------------
Deskos, G., S. Laizet, and M. D. Piggott. “Turbulence-resolving simulations of
wind turbine wakes”. Renewable Energy 134 (2019), pp. 989 –1002. doi: 10.1016/
j.renene.2018.11.084.

Deskos G., S. Laizet, and R. Palacios. “WInc3D: A novel framework for turbulence-resolving 
simulations of wind farm wake interactions”. Wind Energy 23.3 (2020), pp. 779–794.doi:10.1002/we.2458.

### Support or Contact
Georgios.Deskos@nrel.gov 
s.laizet@imperial.ac.uk
r.palacios@imperial.ac.uk

