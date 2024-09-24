Here we have provided a sample simulation. After building lammps as `lmp`, the simulation can be replicated using:

`lmp -in input.in`

This file takes the polymer defined in `init.dat` as the initial file. After runung the code, the output containes `mid_config.dat` which is a relaxed polymer file, `final_config.dat` which contains the final configuration for the polymer, `run_traj` which is the trajectory of the evolution of the polymer and the log files. 

In the default parameters the simulation shows effective acetylation of a polymer strand for 200000 simulation steps. The subsequent analysis on the trajectory file has been shown in `test_analysis.ipynb` which shows the evolution of the polymer in time.

