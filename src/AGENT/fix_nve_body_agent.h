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
  double dtq;                    // timestep length
  double growth_rate;            // expectation of growth rate, unit is 1/[T]
  double growth_standard_dev;    // standard derivation of growth rate
  double L_max;                  // maximum length for proliferation
  double nu_0;                   // viscous coefficient of ambient environments
  double kLH;                    // probability from low cdGMP to high cdGMP
  double kHL;                    // probability from high cdGMP to low cdGMP
  double noise_level;            // pre-defined noise level applying on both force and moment vector
  double coeff_nu_0_xy;          // fold of 3D env viscosity difference on x-y direction
  double coeff_nu_0_z;           // z_direction fold of env viscosity difference
  double z_damp_height;          // apply nu_0 difference above this height

  class AtomVecBody *avec;

  void grow_single_body(int, double);                                         // grow a single cell in a given timestep
  void proliferate_single_body(int);                                          // check length and proliferate a cell             
  void body2space(double*, double*, double*);                                 // convert a vector in body frame to space frame
  void translate_single_body(int, double*);                                   // give a displacement to a cell
  void apply_damping_force(int, double*, double**, double**);                 // environmental viscous force
  double radius(double*, int);                                                // return the radius of a cell
  double length(double*);                                                     // return the length of a cell
  void set_force(int, double, double, double, double, double, double);        // manually set the force and the torque for a given cell
  void add_noise(double*, double*, double);                                   // add random noise to force and moment vectors
  void copy_atom(int, int);                                                   // copy atom information from local index i to j

  void read_params(int, char **);                                                         // read parameters from input script
//   void find_maxid();
};

}    // namespace LAMMPS_NS
#endif
#endif
