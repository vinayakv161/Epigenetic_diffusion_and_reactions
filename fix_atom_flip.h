/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(atom/flip,FixAtomFlip);
// clang-format on
#else

#ifndef LMP_FIX_ATOM_FLIP_H
#define LMP_FIX_ATOM_FLIP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAtomFlip : public Fix {
 public:
  FixAtomFlip(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void pre_exchange();

 private:
  int nevery, seed;
  int ncycles;
  int swap_from, swap_to;


  double nswap_attempts;
  double nswap_successes;
  double probib;


  double energy_stored;

  class RanPark *random_equal;

  class Compute *c_pe;

  int attempt_swap();
  double energy_full();
  int pick_i_swap_atom();
  //int pick_i_swap_atom();
  //int pick_j_swap_atom(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
