/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/body/agent,FixNVEBodyAgent);
// clang-format on
#else

#ifndef LMP_FIX_NVE_BODY_AGENT_H
#define LMP_FIX_NVE_BODY_AGENT_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEBodyAgent : public FixNVE {
 public:
  FixNVEBodyAgent(class LAMMPS *, int, char **);
  void init() override;
  int setmask() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void pre_exchange() override;

 private:
  double dtq;
  class AtomVecBody *avec;

  void grow_single_body(int, double);
  void proliferate_single_body(int);
  void body2space(double*, double*, double*);
  void translate_single_body(int, double*);
  void apply_damping_force(int, double*, double**, double**);
  double radius(double*, int);
  double length(double*);
  void set_force(int, double, double, double, double, double, double);
  void add_noise(double*, double*, double);
  void copy_atom(int, int);
//   void find_maxid();
};

}    // namespace LAMMPS_NS
#endif
#endif
