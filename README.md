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
module load intelmpi (depending on the server)

mpicc -qopenmp -O -I. fdmam-distrib.c nrutil.c -o (exe_file_name)

### RUNING ###
(exe_file_name)

### REQUIREMENT ###

