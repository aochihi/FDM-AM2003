# FDM-AM2003
Finite difference code in C for seismic wave propagation (version Aochi &amp; Madariaga, BSSA, 2003).

This is a finite difference code for seismic wave propagaiton in 3D semi-infinite elastic medium to simulate the ground motions from an earthquake source. 
The formulation is based on the staggard grid wiht 4th order in space and 2nd in time (See the references in Aochi &amp; Madariaga, 2003). 
The program is written in C language, using some libraries from Numerical Recipes. 
The MPI parallelization is built for one-horizontal direcion (X), and OpenMP can be used in each sub-domains. 
A recent version written in C++ is pulished as Ondes3D, availbale from https://bitbucket.org/fdupros/ondes3d/src/master/. 

For any utilization, please refer and cite : 
Aochi, H. and R. Madariaga (2003), The 1999 Izmit, Turkey, earthquake: Non-planar fault structure, dynamic rupture process and strong ground motion, Bulletin of Seismological Society of America, 93, 1249-1266. doi:10.1785/0120020167.

### COMPILING ### 
> module load intelmpi (depending on the server)

> mpicc -qopenmp -O -I. fdmam-distrib.c nrutil.c -o (exe_file_name)

### RUNING ###
> (exe_file_name)

Input files are defined in parameter file. Outputs are written under a directory defined in parameter file. 

### REQUIREMENT ###
>  fdmam-distrib.c (main program)
  
>  "test.prm" (main parameter file, defined as PRM in the main program). To edit. 

>  test.src (souce location file, defined in "test.prm"). To edit. 
  
>  test.hist (source history file, defined in "test.prm"). To edit.
  
>  test.sta (station location file, defined in "test.prm"). To edit. 
  
>  nrutil.c, nrutil.h, nr.h (common files after Numerical Recipes). Don't touch.
  
>  myfdm5-3.h (common file). Don't touch.
  
### Model setting ###
The code calculate the wave propagation in 3D infintie elastic medium (X:EW, Y:NS, Z:Up) with the free surface at Z = 0 m. For simplicity, the 1D layered model can be given in parameter file (test.prm). The numerical parameters should be checked prior to the simulations (no automatic check during run). The earthquake source can be given for any, multiple positions (test.src) wiht complex slip history (test.hist). The outpus are kept for the predefined reveivers positions (test.sta). 

The model parameters are in "test.prm"  

> 4           # 4th order in space. 4 or 2.
>
> -375  625   # X-along grid number XMIN, XMAX
>
> -300  300   # Y-along grid number YMIN, YMAX
>
> -200        # Z-along grid number ZMIN, ZMAX=0 (implicit)
>
> 6000        # total time steps

> ./test1/    # output directory

> test1.src   # source position file

> test1.hist  # source slip history file

> test1.sta   # station position file

> 200   0.01  # grid size (ds) in m, time step (dt) in s 

> 4           # number of 1D layer model

> 0   3220  1780  2330  300 # top depth of layer 1 (km), Vp (m/s), Vs (m/s), rho (kg/m3), Q

> -1  4620  2650  2550  300 # layer 2

> -5  5940  3420  2710  300 # layer 3

> -12 6170  3550  2750  300 # layer 4 (the end)

