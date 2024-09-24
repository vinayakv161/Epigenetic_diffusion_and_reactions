# Epigenetic_diffusion_and_reactions

This is the repository for the simulation code of the manuscript: "Polymer Model Integrates Super-Resolution Imaging and Epigenomic Sequencing to Elucidate the Role of Epigenetic Reactions in Shaping 4D Chromatin Organization." We have built two custom fixes for [`LAMMPS`](https://www.lammps.org/) to execute epigenetic reaction function (fix atom/swap) and diffusion of epigenetic marks (fix kawasaki). 

## How to use:

### Build LAMMPS:

1. Put `fix_kawasaki.cpp/fix_kawasaki.h` and `fix_atom_flip.cpp/fix_atom_flip.h` in `src/MC` folder of `LAMMPS`.
2. Build `LAMMPS` with the `MC` package

### Use in input file:

1. Using fix kawasaki:
`fix ID group-ID kawasaki Nevery Nparticles seed T`

* `ID,group-ID`: could be found in the [fix](https://docs.lammps.org/fix.html) documentation of `LAMMPS`.
* `Nstep`: attempt kawasaki every this many steps
* `Nparticles`: number of group atims to consider for kawasaki
* `seed`: random # seed (positive integer)
* `T`: Temperature of the kawasaki metropolis criteria

2. Using fix atom/flip:
`fix ID group-ID atom/flip Nevery Nparticles seed type1 type2 prob`

* `ID, group-ID`: could be found in the [fix](https://docs.lammps.org/fix.html) documentation of `LAMMPS`.
* `Nevery`: attempt atom/flip every this many steps
* `Nparticles`: number of group atims to consider for kawasaki
* `seed`: random # seed (positive integer)
* `type1, type2, prob`: The transition probability of type1 atom to type2 atoms is `prob`

### Example

