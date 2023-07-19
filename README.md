## Code for solving spatially heterogeneous coagulation problem sets in 1D with the source of monomers.

The monomer source at the left boundary allows us to study the stationary solution

* Parameters:
`h` - particle size step, 
`dt` - time step, 
`J` - strength of the source, 
`TIME` - total number of iterations, 
`MOD` - data dump periodicity, 
`max_particle_size` - total number of particle sizes, 
`max_x` - total number of space steps, 
`dx` - space step, 
`vel_coefs` - array of advection coefficients, 
`dif_coefs` - array of diffusion coefficients, 
`initial_layer` - array with the initial condition (default is zero).

* Monomer source is presented as a boundary condition in the file __main.cpp__ (line **290** - `rhs`). Also see lines **338-345** (this shows the values of the first and the last rows of tridiagonal matrix. Using it we can set different boundary condtition types, e.g. Neumann on the left: `[0+m*max_x]` and Dirichlet on the right: `[(max_x-1)+m*max_x]`.

* Data is dumped every `MOD` iterations in time. See file __baikal_plots.ipynb__. It shows how the data is read using Python.

* The repository includes 5 configuration files (__configs__ folder):
0. model on July 23rd with smooth ozone source
1. model on July 23rd with direct ozone source
2. model on July 25th with smooth ozone source
3. model on July 25rd with direct ozone source
4. countinuous 72 hour model from July 23rd to July 26th with direct ozone source

To launch the necessary simulation pass the respective configuration file as the command line argument to the executive file.

* Command line arguments. For instance, if one launches the executable with `solver.exe output.txt 0 I.txt configs/input_1.txt`, then __output.txt__ file will apprear with the aerosol data and __I.txt__ with the ozone data. The configurations will be taken from __configs/input_1.txt__.
