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
   Contributing authors: Vinayak (UPenn, ShenoyLab)
------------------------------------------------------------------------- */

#include "fix_kawasaki.h"

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
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "random_park.h"
#include "region.h"
#include "update.h"
#include <iostream>

#include <cmath>
#include <cctype>
#include <cfloat>
#include <cstring>
#include <tuple>

using namespace LAMMPS_NS;
using namespace FixConst;

FixKawasaki::FixKawasaki(LAMMPS *lmp, int narg, char **arg):
  Fix(lmp, narg, arg), list(nullptr),
  random_equal(nullptr),  c_pe(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal fix kawasaki command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  ncycles = utils::inumeric(FLERR,arg[4],false,lmp);
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
  double temperature = utils::numeric(FLERR,arg[6],false,lmp);

  if (nevery <= 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (ncycles < 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (seed <= 0) error->all(FLERR,"Illegal fix atom/swap command");
  if (temperature <= 0.0) error->all(FLERR,"Illegal fix atom/swap command");

  //computing inverse temperature
  beta = 1.0/(force->boltz*temperature);

  //random number generator
  random_equal = new RanPark(lmp,seed);

  //set up reneighboring
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // zero out counters
  nswap_attempts = 0.0;
  nswap_successes = 0.0;
  

}

int FixKawasaki::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

void FixKawasaki::init()
{
  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];
  //std::cout << "init" << "\n";
  int *type = atom->type;

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half=0;
  neighbor->requests[irequest]->full=1;
  neighbor->requests[irequest]->occasional = 1;
  //neighbor->requests[irequest]->cut=1;
  //neighbor->requests[irequest]->cutoff=1.1+neighbor->skin;
}


void FixKawasaki::init_list(int /*id*/, NeighList *ptr)
{
  //std::cout << "-1" << "\n";
  list = ptr;
}

void FixKawasaki::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);

  energy_stored = energy_full();
  //std::cout << "0" << "\n";
  int nsuccess = 0;
  for (int k=0;k<ncycles;k++) nsuccess += attempt_swap();

  nswap_attempts += ncycles;
  nswap_successes += nsuccess;

  next_reneighbor = update->ntimestep + nevery;
}

int FixKawasaki::attempt_swap()
{
   //std::cout << "1" << "\n";
   double energy_before = energy_stored;
   //std::cout << "En before calculated" << "\n";
   int i,j;
   pick_i_j_swap_atom(i,j);
   //std::cout << "picked i and j" << "\n";
   int itype = atom->type[i];
   int jtype = atom->type[j];
   
   if (itype==3){
     return 0;
   }

   if (jtype==3){
     return 0;
   }

   if (i==-1){
      if(j==-1) return 0;
   } 

   if (itype==jtype){
      return 0;
   }
   
   atom->type[i] = jtype;
   atom->type[j] = itype;
   //std::cout << "swapped i and j" << "\n";
   double energy_after = energy_full();
   //std::cout << "energy after" << "\n";
   if (random_equal->uniform() <
      exp(beta*(energy_before - energy_after))) {
         //std::cout << "YES Atom list" << atom->nlocal << "\n";
    energy_stored = energy_after;
    return 1;
   }
   //std::cout << "swapping complete" << "\n";
   // swap not accepted, return 0
   // restore the swapped itype & jtype atoms
   // do not need to re-call comm->borders() and rebuild neighbor list
   //   since will be done on next cycle or in Verlet when this fix finishes

   
   atom->type[i] = itype;
   atom->type[j] = jtype;

   return 0;
}

double FixKawasaki::energy_full()
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


/*int FixKawasaki::pick_i_swap_atom()
{
  int nlocal = atom->nlocal;
  int i=-1;
  i = static_cast<int>(nlocal*random_equal->uniform());

  return i;
  
}

int FixKawasaki::pick_j_swap_atom(int i)
{
  int ij;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int nlocal = atom->nlocal;

  int j = -1;
  
  neighbor->build_one(list,1);
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  int gotj=1;
  int counter=0;
  jlist=firstneigh[i];
  int jnum=numneigh[i];
  
  int j_index = static_cast<int>(jnum*random_equal->uniform());

  j = jlist[j_index];

  return j;
  
}*/

void FixKawasaki::pick_i_j_swap_atom(int &a, int &b)
{
  int nlocal = atom->nlocal;
  int i=-1;
  int j=-1;
  double rand = random_equal->uniform();
  int randindex = static_cast<int>(nlocal*rand);

  int inum,*ilist,*jlist,*numneigh,**firstneigh;
  
  neighbor->build_one(list,1);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  i = ilist[randindex];
  jlist=firstneigh[i];
  int jnum=numneigh[i];

  if (jnum==0){
     a=-1;
     b=-1;
     return;
  }
  
  int j_index = static_cast<int>(jnum*random_equal->uniform());

  j = jlist[j_index];

  a = i;
  b = j;
  
}
