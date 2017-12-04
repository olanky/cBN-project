# cBN-project

## Description

This set of code is used to simulate the bombardment of particles (NH molecules, N or H atoms) toward the center of mass (COM) of the target cluster, working with DFTB+. 

format.cpp code will read in the coordinates of previous bombardment. And perform following steps:

1. Remove any particle leaving cluster.
2. Translate the cluster so that the COM of cluster is at the origin (0,0,0).
3. Deposite a new particle at the sperical surface with radius of 3 nm from COM of target cluster, directing new particle toward the target cluster, with thermal velocity corresponding to simulation temperature.

The bash script is used to achieve continous bombardment with the new particles.

## Usage

1. Compile format.cpp code:
> $ icpc format.cpp -o format
2. Set up parameters consistently in DFTB input and bash script, including time period for each bombardment and simulation temperature.
3. Specify type of projectile particles in bash script.

