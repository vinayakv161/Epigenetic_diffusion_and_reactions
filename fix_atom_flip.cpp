// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Vinayak (Snenoy Lab, UPenn)
------------------------------------------------------------------------- */

#include "fix_atom_flip.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"

#include <cmath>
#include <cctype>
#include <cfloat>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAtomFlip::FixAtomFlip(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random_equal(nullptr), c_pe(nullptr)
{
  //if (narg < 20) error->all(FLERR,"Illegal fix atom/swap command");

  // required args

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  ncycles = utils::inumeric(FLERR,arg[4],false,lmp);
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
  swap_from = utils::inumeric(FLERR,arg[6],false,lmp);
  swap_to = utils::inumeric(FLERR,arg[7],false,lmp);
  double prob = utils::numeric(FLERR,arg[8],false,lmp);
  //double temperature = utils::numeric(FLERR,arg[6],false,lmp);

  probib = prob;

  if (nevery <= 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (ncycles < 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (seed <= 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (probib <= 0.0) error->all(FLERR,"Illegal fix atom/swap command");

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters

  nswap_attempts = 0.0;
  nswap_successes = 0.0;

}

/* ---------------------------------------------------------------------- */

int FixAtomFlip::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAtomFlip::init()
{
  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  int *type = atom->type;
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo swaps
------------------------------------------------------------------------- */

void FixAtomFlip::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // insure current system is ready to compute energy

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);

  // energy_stored = energy of current state
  // will be updated after accepted swaps

  energy_stored = energy_full();

  // attempt Ncycle atom swaps

  int nsuccess = 0;
  for (int i = 0; i < ncycles; i++) nsuccess += attempt_swap();

  // udpate MC stats

  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  next_reneighbor = update->ntimestep + nevery;
}

/* ----------------------------------------------------------------------
   attempt a semd-grand swap of a single atom
   compare before/after energy and accept/reject the swap
   NOTE: atom charges are assumed equal and so are not updated
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   attempt a swap of a pair of atoms
   compare before/after energy and accept/reject the swap
------------------------------------------------------------------------- */

int FixAtomFlip::attempt_swap()
{
  
  // pre-swap energy

  double energy_before = energy_stored;

  // pick a random pair of atoms
  // swap their properties

  int i = pick_i_swap_atom();

  if (atom->type[i]==3){
    return 0;
  }
  
  if (atom->type[i]>=0){
    atom->type[i] = swap_to;
  }
  // if unequal_cutoffs, call comm->borders() and rebuild neighbor list
  // else communicate ghost atoms
  // call to comm->exchange() is a no-op but clears ghost atoms

  // post-swap energy

  double energy_after = energy_full();

  // swap accepted, return 1
  // if ke_flag, rescale atom velocities

  if (random_equal->uniform() <= probib) {
    energy_stored = energy_after;
    return 1;
  }

  // swap not accepted, return 0
  // restore the swapped itype & jtype atoms
  // do not need to re-call comm->borders() and rebuild neighbor list
  //   since will be done on next cycle or in Verlet when this fix finishes

  if (i >= 0) {
    atom->type[i] = swap_from;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixAtomFlip::energy_full()
{
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  if (modify->n_post_force) modify->post_force(vflag);

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixAtomFlip::pick_i_swap_atom()
{
  int i = -1;
  int nlocal = atom->nlocal;
  i = static_cast<int>(nlocal*random_equal->uniform());

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */
