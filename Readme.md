# Perturbative gradient flow coupling

Perturbative_gf_coupling is a C++ package which computes the gradient flow coupling normalisation at leading order in perturbation theory (O(g_0^2)).

In finite volume, with Schroedinger Functional boundary conditions, a non perturbative definition of the strong coupling constant requires a normalisation factor. Depending of the particular observable chosen, it can be computed in perturbation theory. Here we consider the GF coupling defined in terms of the so-called action density E(t,x) (see reference by M. Luescher [arXiv:1404.5930]).
The first numerical result were presented at the lattice conference 2016 (see reference by A. R. and S. Sint [arXiv:1612.07047]), further studies are presented in my PhD thesis (to appear).
The improvement coefficient for the initial condition in the flow equation is discussed in the reference by A. Ramos and S. Sint [arXiv:1508.05552].


## Installation
### Dependences
In order to run the code, you need to have a g++ compiler and to download the armadillo library (http://arma.sourceforge.net/download.html).

### Build
Use the Makefile to generate the executable

```bash
cd build
make .
```

## Usage

```
Nmag ob fl ac c L/a T/a cb tau
EXAMPLE: Nmag plcl z pl 0.3 16 16 0.1 0
ob = observable, options: pl,cl,plcl,lw
fl=flow, options: fl = w,z,s
ac=action, options: ac = pl,lw
c= sqrt(8t)/L  
L/a = spatial lattice size
T/a = temporal lattice size
cb = imporvement coeff. in the initial condition for the flow
tau =(shift in flow time)
```

```
Nelec eps ob fl ac c L/a T/a c_b
EXAMPLE: Nelec 0.01 imp z pl 0.3 8 8 0.02
eps = parameter to compute numerical derivative
c= sqrt(8t)/L  
ob = observable, options: pl (plaq not symm),pls (plaq symm), cl (clover (symm by def)),plcl ()
fl=flow, options: fl = w,z,s
ac=action, options: ac = pl,lw
L/a = spatial lattice size
T/a = temporal lattice size
c_b = imporvement coeff. in the initial condition for the flow
```
You may want to redirect the output to collect this information in the folder data.

## Contact
For any problems, suggestions, ideas, please contact us: distinguibile@gmail.com 

## License
GNU GENERAL PUBLIC LICENSE
