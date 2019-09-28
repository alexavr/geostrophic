*Detailed description is coming soon...*

# Description

Computes *geostrophic* velocity using geostrophic balance. Additionally computes *ageostrophic* velocity as a vector subtraction of total velocity and geosprophic velocity. 

Canonically the total velocity is the decomposition of geostrophic and ageostrophic velocity:

$\vec{V}=\vec{V_g}+\vec{V_a}$

Where the geostrophic velocity corresponds to synoptic motions and can be trivially derived using geostrophic balance

$\vec{V_g}=\frac{1}{\rho}\Delta{P}$

In order to gain the ageostrophic velocity as the mesoscale marker we compute geostrophic velocity first and then subtract it from total velocity.

For WRF output there is an option to compute the geostrophic velocity on model levels using fluctuation fields instead of total variables using numerical schemes taken from WRF. This increases the accuracy. You need to run this code upon true WRF output.

## Requirements

[NetCDF](https://www.unidata.ucar.edu/).

## Building 

Copy to your machine:

```bash
git clone https://github.com/alexavr/geostrophic.git
cd geostrophic
vi Makefile
make
```

You might need to edit `Makefile` :

```bash
NETCDF=/opt/netcdf4-hdf5       # path to your NetCDF
FC=ifort                       # your FORTRAN compiler
```

Cleaning in case anything went wrong:

```bash
make clean
```

## Data preprocessing

If you're planing to compute geostrophic upon WRF output than you don't have to prepare anything.

If you're planing to compute upon other data (say reanalyse) only simple (or no-fluctuation) method is available. The source file has to contain *u,v,p,rho* variables along with *lat*, *lon* (or *latitude*, *longitude*) fields. This preparation could be accomplished by using the [cdo](https://code.mpimet.mpg.de/projects/cdo/) library. 

No moving grids supported yet.

## Configuration

The configuration is done using the file *namelist.geo*.



The domain could be regional (`regional = .TRUE.`) or global (periodic boundary, `regional = .FALSE.`). You need also set the time interval using `stime` and `etime` (this interval has to be present in the source file!).

The particles could be initiated manually (`pt_grid = .FALSE.`) in ASCII file named *init_particles.dat* (see *init_particles_Reykjavik.dat* as an example). Or evenly spaced allover the domain (`pt_grid = .TRUE.`), `pt_step ` is number of steps between particles in both directions, `pt_height` is elevation in meters. 

Accuracy of the computations. There are two schemes provided: 

1. Coarse scheme (`accuracy = .FALSE.`). Predicts particle location using time step from the source file. Fast but very inaccurate. Recommended only for test purposes. 
2. Fine scheme (`accuracy = .TRUE.`). Interpolates data between source file time step. This option controlled by `timestep` (in minutes).

## Running

```bash
./geostrophic.exe src_file.nc
```

where *src_file.nc* comes from the "Data reprocessing".



