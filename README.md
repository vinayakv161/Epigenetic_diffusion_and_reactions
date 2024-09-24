# Epigenetic_diffusion_and_reactions

This is the repository for the simulation code of the manuscript: "Polymer Model Integrates Super-Resolution Imaging and Epigenomic Sequencing to Elucidate the Role of Epigenetic Reactions in Shaping 4D Chromatin Organization." We have built two custom fixes for 'LAMMPS' to execute epigenetic reaction function (fix atom/swap) and diffusion of epigenetic marks (fix kawasaki). 

## How to use:

### Build LAMMPS:

1. Put fix_kawasaki.cpp/fix_kawasaki.h and fix_atom_flip.cpp/fix_atom_flip.h in src/MC folder of LAMMPS.
2. Build LAMMPS with the MC package

### Use in input file:

1. Using fix kawasaki:
fix ID group-ID kawasaki

3. Using fix atom/flip:
