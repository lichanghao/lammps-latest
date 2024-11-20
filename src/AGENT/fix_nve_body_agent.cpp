// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Agent-based simulation for bacteria biofilms
   Author: Changhao Li (changhaoli1997@gmail.com)
   Last updated: 11/17/2024
------------------------------------------------------------------------- */

#include <cmath>
#include "fix_nve_body_agent.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "error.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INITIAL_MASS 1e-3
#define INITIAL_INERTIA_x 5e-4
#define INITIAL_INERTIA_y 1.2708e-4
#define INITIAL_INERTIA_z 1.2708e-4
#define RANDOM_SEED 2189634
#define FORCE_RENEIGHBOR_INTERVAL 1

// #define FIX_NVE_BODY_AGENT_DEBUG
// #define DEBUG_INTERVAL 2000

/* ---------------------------------------------------------------------- */

FixNVEBodyAgent::FixNVEBodyAgent(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) 
{ 
  // set fix parent class tags
  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  time_integrate = 1;

  // read parameters from input files
  read_params(narg, arg);

  // random generator (seed = RANDOM_SEED + processor_id)
  random = new RanPark(lmp, RANDOM_SEED + comm->me);

  // initiate peratom vector for growth rates, Gaussian distribution ~ N(growth_rate, growth_standard_dev)
  nmax = atom->nmax;
  grow_arrays(nmax);
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    // if (mask[i] & groupbit)
      growth_rates_all[i] = random->gaussian() * growth_standard_dev + growth_rate;
  }
  atom->add_callback(Atom::GROW);

  // initiate the image flag for all atoms as 0, because somehow the original body package did not do it
  for (int i = 0; i < nlocal; i++) atom->image[i] = 0;

  // find maximum id across all processors
  find_maxid();
}

/* ---------------------------------------------------------------------- */

void FixNVEBodyAgent::init()
{
  avec = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
  if (!avec) error->all(FLERR,"Fix nve/body/agent requires atom style body");

  avec_hybrid = dynamic_cast<AtomVec *>(atom->style_match("hybrid"));
  if (!avec) avec_hybrid = avec;

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;

  // check that all particles are bodies
  // no point particles allowed

  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (body[i] < 0) error->one(FLERR,"Fix nve/body/agent requires bodies");

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

FixNVEBodyAgent::~FixNVEBodyAgent()
{
  delete random;
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(growth_rates_all);
}

/* ---------------------------------------------------------------------- */

int FixNVEBodyAgent::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_EXCHANGE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- 
  do the first half of verlet integration, apply damping and noise
---------------------------------------------------------------------- */

void FixNVEBodyAgent::initial_integrate(int /*vflag*/)
{
  double dtfm;
  double omega[3];
  double *quat,*inertia;

  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];

      // apply 3D ambient damping forces (F = -mu * v, T = -mu * omega)
      inertia = bonus[body[i]].inertia;
      quat = bonus[body[i]].quat;
      MathExtra::mq_to_omega(angmom[i], quat, inertia, omega);
      apply_damping_force(i, omega, f, torque);
      
      // apply infinitesimal noise to break symmetry
      add_noise(f[i], torque[i], noise_level);

      // update velocity by full step, displacement by 1/2 step
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step
      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion
      MathExtra::mq_to_omega(angmom[i], quat, inertia, omega);
      MathExtra::richardson(quat, angmom[i], omega, inertia, dtq);
    }

    determine_next_reneighbor();
}


/* ---------------------------------------------------------------------- 
  insert new particles
---------------------------------------------------------------------- */

void FixNVEBodyAgent::pre_exchange()
{ 
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  // just return if should not be called on this timestep
  if (next_reneighbor != update->ntimestep) return;

  // clear atom maps and ghost atom information
  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // loop over all cells and determine which one is dividing
  int nadded = 0;
  for (int i = 0; i < nlocal; i++) {
    bool any_division = false;
    if (mask[i] & groupbit) {
      proliferate_single_body(i, any_division);
      if (any_division) {
        nadded++; 
        // note: currently tag is randomly initialized from INT_MAX to the initial maxtag_all
        // tag > 1e6 will cause segmentation fault if you use array style atom map
        // should use hash style atom map instead to avoid this problem
        tagint newtag = maxtag_all + static_cast<tagint>(random->uniform() * (MAXTAGINT - maxtag_all));
        atom->tag[atom->nlocal - nadded] = newtag;
      }
    }
  }

  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init(1);
    atom->map_set();
  }
}


/* ---------------------------------------------------------------------- 
  do the second half verlet integration, grow all rod-like bodies
---------------------------------------------------------------------- */

void FixNVEBodyAgent::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      grow_single_body(i, growth_rates_all[i]);
    }
}


/* ----------------------------------------------------------------------
  see if this body needs proliferation
  if dividing, insert a new daughter body
---------------------------------------------------------------------- */

void FixNVEBodyAgent::proliferate_single_body(int ibody, bool &is_dividing)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  double *rmass = atom->rmass;
  double **angmom = atom->angmom;
  int *body = atom->body;
  double **x = atom->x;

  double *p = bonus[body[ibody]].dvalue;
  double L = length(p);                             // length of the rod
  double r = radius(p, 2);                          // rounded radius of the rod

  if (L >= L_max) {
    is_dividing = true;
    // calculate new center locations of two daughter bodies
    // c1 is the relative displacements of the first daughter body, c2 is the second
    #ifdef FIX_NVE_BODY_AGENT_DEBUG
    printf("proliferating the cell %d with rounded radius %f, current timestep: %d\n", ibody ,bonus[body[ibody]].dvalue[6+2], update->ntimestep);
    #endif
    double temp[6];
    for (int j = 0; j < 6; j++)
      temp[j] = bonus[body[ibody]].dvalue[j] * (L / 2 + r) / L;
    double cc[6];
    body2space(temp, bonus[body[ibody]].quat, cc);
    double *c1 = cc;
    double *c2 = cc + 3;

    double x_old[3];
    for (int j = 0; j < 3; j++) {
      x_old[j] = x[ibody][j];
    }

    // mother body -> the first daughter body
    double prolif_ratio = (L / 2 - r) / L;
    for (int j = 0; j < 6; j++)
      bonus[body[ibody]].dvalue[j] *= prolif_ratio;
    bonus[body[ibody]].dvalue[6+2] *= prolif_ratio;
    translate_single_body(ibody, c1);                  // translate the mother body

    // insert the second daughter body
    // avec->add_body(ibody);
    for (int j = 0; j < 3; j++)
    {
      // x[new_body_index][j] = x_old[j] + c2[j];      // set up the center location of the second daughter body
      c2[j] = c2[j] + x_old[j];
    }
    atom->avec->create_atom(atom->type[ibody], c2);
    int new_body_index = atom->nlocal - 1;             // append the new cell to the last
    copy_atom(ibody, new_body_index);                  // copy the atom information, such as velocity and angular momentum
    body[new_body_index] = 0;
    int int_temp[3] = {2, 0, 0};
    double double_temp[13] = {INITIAL_INERTIA_x, INITIAL_INERTIA_y, INITIAL_INERTIA_z, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              -1.250000, 0.000000, 0.000000, 1.250000, 0.000000, 0.000000, 2.000000};
    avec->data_body(new_body_index, 3, 13, int_temp, double_temp);
    avec->deep_copy_bonus(body[ibody], body[new_body_index]);

    // forcing no net external forces
    set_force(ibody, 0, 0, 0, 0, 0, 0);
    set_force(new_body_index, 0, 0, 0, 0, 0, 0);

    // reset the mass and rotational inertia
    rmass[new_body_index] = INITIAL_MASS;
    rmass[ibody] = INITIAL_MASS;
    bonus[body[ibody]].inertia[0] = INITIAL_INERTIA_x;
    bonus[body[ibody]].inertia[1] = INITIAL_INERTIA_y;
    bonus[body[ibody]].inertia[2] = INITIAL_INERTIA_z;
    bonus[body[new_body_index]].inertia[0] = INITIAL_INERTIA_x;
    bonus[body[new_body_index]].inertia[1] = INITIAL_INERTIA_y;
    bonus[body[new_body_index]].inertia[2] = INITIAL_INERTIA_z;

    // reallocate the per-atom vector if nmax is changed
    if (nmax < atom->nmax) {
      nmax = atom->nmax;
      grow_arrays(nmax);
    }
    growth_rates_all[new_body_index] = random->gaussian() * growth_standard_dev + growth_rate;
  }
}


/* ----------------------------------------------------------------------
  grow single body, index is ibody
  growth_ratio is alpha in 2018 Nat. Phys. Paper (Farzan Beroz et al.)
---------------------------------------------------------------------- */

void FixNVEBodyAgent::grow_single_body(int ibody, double growth_rate)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  double *rmass = atom->rmass;
  int *body = atom->body;
  double r = radius(bonus[body[ibody]].dvalue, 2);
  double L = length(bonus[body[ibody]].dvalue);

  if (bonus[body[ibody]].ndouble != 10) {printf("warning: ndouble != 10, here ndouble = %d\n", bonus[body[ibody]].ndouble);}

  // add Gaussian random noise to growth rate: current sigma = 0.2<alpha>

  // dL/L = alpha * (4/3*R + L)/L * dt
  // dV/V = alpha * dt
  double length_ratio = 1 + 2 * dtf * growth_rate * (4.0 / 3.0 * r + L) / L;
  double growth_ratio = 1 + 2 * dtf * growth_rate;

  // update the coords of vertices
  for (int j = 0; j < 6; j++) {
    bonus[body[ibody]].dvalue[j] *= length_ratio; // coords of vertices (in body frame)
  }
  bonus[body[ibody]].dvalue[6+2] *= length_ratio; // enclosing radius (not rounded radius)

  // mass and rotational inertia are currently neglected, because the model is running under overdamped settings
  rmass[ibody] *= growth_ratio;                   // mass ~ V
  for (int j = 0; j < 3; j++) {
    bonus[body[ibody]].inertia[j] *= std::pow(growth_ratio, 3); // rotational inertia ~ m*L^2
  }
}


/* ----------------------------------------------------------------------
  transfer the body-frame relative displacement to space-frame
  specifially for rods
  don't forget to add the coords of the rod center!
---------------------------------------------------------------------- */

void FixNVEBodyAgent::body2space(double *coords, double *quat, double *space_coords)
{
  double p[3][3];
  // translate quaterions to rotation matrix R
  MathExtra::quat_to_mat(quat, p);
  for (int m = 0; m < 2; m++)
  {
    // d_space = R * d_body, do it for each vertex
    MathExtra::matvec(p, &coords[3 * m], &space_coords[3 * m]);
  }
}


/* ----------------------------------------------------------------------
  translate single body
---------------------------------------------------------------------- */

void FixNVEBodyAgent::translate_single_body(int ibody, double *vec)
{
  double **x = atom->x;
  for (int j = 0; j < 3; j++)
    x[ibody][j] += vec[j];
}


/* ----------------------------------------------------------------------
  apply 3D viscous drags, including translational and rotational parts
---------------------------------------------------------------------- */

void FixNVEBodyAgent::apply_damping_force(int ibody, double *omega, double **f, double **torque)
{
  AtomVecBody::Bonus *bonus = avec->bonus;
  double *v = (atom->v)[ibody];
  int *body = atom->body;
  double L = length(bonus[body[ibody]].dvalue);
  double R = radius(bonus[body[ibody]].dvalue, 2);

  double temp_nu_0 = nu_0;

  // adding damping force, applying on mass center
  f[ibody][0] += -temp_nu_0 * (L+4.0/3.0*R) * v[0];
  f[ibody][1] += -temp_nu_0 * (L+4.0/3.0*R) * v[1];
  f[ibody][2] += -temp_nu_0 * (L+4.0/3.0*R) * v[2];

  // adding damping moment
  torque[ibody][0] += -1.0 / 6.0 * temp_nu_0 * omega[0] * std::pow(L+4.0/3.0*R, 3);
  torque[ibody][1] += -1.0 / 6.0 * temp_nu_0 * omega[1] * std::pow(L+4.0/3.0*R, 3);
  torque[ibody][2] += -1.0 / 6.0 * temp_nu_0 * omega[2] * std::pow(L+4.0/3.0*R, 3);
  
  // debug code
  #ifdef FIX_NVE_BODY_AGENT_DEBUG
  if (update->ntimestep % DEBUG_INTERVAL == 0) {
    printf("apply_damping_force for body %d\n", ibody);
    printf("damping moment: %f %f %f\n", torque[ibody][0], torque[ibody][1], torque[ibody][2]);
    printf("damping force: %f %f %f\n", f[ibody][0], f[ibody][1], f[ibody][2]);
    printf("omega: %f %f %f\n", omega[0], omega[1], omega[2]);
    printf("inertia: %f %f %f\n", bonus[body[ibody]].inertia[0], bonus[body[ibody]].inertia[1], bonus[body[ibody]].inertia[2]);
  }
  #endif
}


/* ----------------------------------------------------------------------
  Adding noise to force vector and moment vector
---------------------------------------------------------------------- */

void FixNVEBodyAgent::add_noise(double *f, double *mom, double given_noise_level)
{
  f[0] += 1 * given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  f[1] += 1 * given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  f[2] += 1 * given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[0] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[1] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
  mom[2] += given_noise_level * (rand() / (static_cast<double>(RAND_MAX)) - 0.5);
}


/* ----------------------------------------------------------------------
  return rounded radius for rod
---------------------------------------------------------------------- */

double FixNVEBodyAgent::radius(double *data, int nvert)
{
  return data[nvert * 3 + 2 + 1];
}


/* ----------------------------------------------------------------------
  return length for rod
---------------------------------------------------------------------- */

double FixNVEBodyAgent::length(double *data)
{
  return sqrt(std::pow(data[3] - data[0], 2) + std::pow(data[4] - data[1], 2) + std::pow(data[5] - data[2], 2));
}


/* ----------------------------------------------------------------------
  manually set the force and torque for the give body
---------------------------------------------------------------------- */

void FixNVEBodyAgent::set_force(int ibody, double fx, double fy, double fz, double tx, double ty, double tz)
{
  double *f = atom->f[ibody];
  double *torque = atom->torque[ibody];

  f[0] = fx;
  f[1] = fy;
  f[2] = fz;
  torque[0] = tx;
  torque[1] = ty;
  torque[2] = tz;
}


/* ----------------------------------------------------------------------
  copy the atom information from ibody to jbody
  location information is not copied
  tag is not copied and initialized as 0
  tag needs to be specially taken care of, especially under parallel settings
---------------------------------------------------------------------- */

void FixNVEBodyAgent::copy_atom(int ibody, int jbody)
{
  // atom->tag[jbody] = 0;
  atom->type[jbody] = atom->type[ibody];
  // atom->x[jbody][0] = atom->x[ibody][0];
  // atom->x[jbody][1] = atom->x[ibody][1];
  // atom->x[jbody][2] = atom->x[ibody][2];
  atom->mask[jbody] = atom->mask[ibody];
  atom->image[jbody] = atom->image[ibody];
  atom->v[jbody][0] = atom->v[ibody][0];
  atom->v[jbody][1] = atom->v[ibody][1];
  atom->v[jbody][2] = atom->v[ibody][2];
  atom->radius[jbody] = atom->radius[ibody];
  atom->rmass[jbody] = atom->rmass[ibody];
  atom->angmom[jbody][0] = atom->angmom[ibody][0];
  atom->angmom[jbody][1] = atom->angmom[ibody][1];
  atom->angmom[jbody][2] = atom->angmom[ibody][2];

  // debug code
  #ifdef FIX_NVE_BODY_AGENT_DEBUG
  printf("copy_atom() from local index %d to %d\n", ibody, jbody);
  printf("location of %d: %f %f %f\n", ibody, atom->x[ibody][0], atom->x[ibody][1], atom->x[ibody][2]);
  printf("location of %d: %f %f %f\n", jbody, atom->x[jbody][0], atom->x[jbody][1], atom->x[jbody][2]);
  printf("image flag: %d\n", atom->image[jbody]);
  printf("atom tag: %d\n", atom->tag[jbody]);
  printf("velocity: %f %f %f\n", atom->v[jbody][0], atom->v[jbody][1], atom->v[jbody][2]);
  printf("angular momentum: %f %f %f\n", atom->angmom[jbody][0], atom->angmom[jbody][1], atom->angmom[jbody][2]);
  #endif
}


/* ----------------------------------------------------------------------
  grow peratom arrays for cell-wise growth rates, etc.
------------------------------------------------------------------------- */

void FixNVEBodyAgent::grow_arrays(int n)
{
  memory->grow(growth_rates_all, n, "fix/nve/body/agent:growth_rates_all");
  vector_atom = growth_rates_all;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixNVEBodyAgent::memory_usage()
{
  double bytes = (double) atom->nmax * 1 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixNVEBodyAgent::copy_arrays(int i, int j, int /*delflag*/)
{
  growth_rates_all[j] = growth_rates_all[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixNVEBodyAgent::set_arrays(int i)
{
  growth_rates_all[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixNVEBodyAgent::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = growth_rates_all[i];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixNVEBodyAgent::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  growth_rates_all[nlocal] = buf[n++];
  return n;
}

/* ----------------------------------------------------------------------
  read LAMMPS fix parameters from input script
------------------------------------------------------------------------- */

void FixNVEBodyAgent::read_params(int narg, char **arg)
{
  // default values
  growth_rate = 0;
  growth_standard_dev = 0;
  L_max = 0;
  nu_0 = 0;
  kLH = 0;
  kHL = 0;
  noise_level = 1e-8;
  coeff_nu_0_xy = 1;
  coeff_nu_0_z = 1;
  z_damp_height = 0;

  if (narg != 3 && narg != 6 && narg != 8 && narg != 10 && narg != 12 && 
      narg != 14 && narg != 16 && narg != 18 && narg != 20 && narg != 22) {
    error->all(FLERR, "Invalid fix nve/body/agent command, incorrect number of input parameters");
  }

  // read values from LAMMPS input file
  for (int i = 0; i < narg - 1; i++)
  {
    if (strcmp(arg[i], "grow") == 0)
    {
      growth_rate = utils::numeric(FLERR, arg[i + 1], false, lmp);
      growth_standard_dev = utils::numeric(FLERR, arg[i + 2], false, lmp);
    }
    if (strcmp(arg[i], "maxlen") == 0)
    {
      L_max = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "damp") == 0)
    {
      nu_0 = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "noise") == 0)
    {
      noise_level = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "kLH") == 0)
    {
      kLH = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "kHL") == 0)
    {
      kHL = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "c_nu_0_xy") == 0)
    {
      coeff_nu_0_xy = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "c_nu_0_z") == 0)
    {
      coeff_nu_0_z = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
    if (strcmp(arg[i], "z_damp_height") == 0)
    {
      z_damp_height = utils::numeric(FLERR, arg[i + 1], false, lmp);
    }
  }

  if (comm->me == 0) {
    printf("\n------------------ Fix_NVE_Body_Agent Parameters ------------------\n");
    printf("growth_rate = %f\n", growth_rate);
    printf("growth_standard_dev = %f\n", growth_standard_dev);
    printf("L_max = %f\n", L_max);
    printf("nu_0 = %f\n", nu_0);
    printf("noise_level (normalized by 1e-8) = %f\n", noise_level/1e-8);
    printf("kLH = %f\n", kLH);
    printf("kHL = %f\n", kHL);
    printf("coeff_nu_0_xy = %f\n", coeff_nu_0_xy);
    printf("coeff_nu_0_z = %f\n", coeff_nu_0_z);
    printf("z_damp_height = %f\n", z_damp_height);
    printf("-------------------------------------------------------------------\n\n");
  }
}

/* ----------------------------------------------------------------------
  loop over all atoms to determine if there is a need to reneighboring
------------------------------------------------------------------------- */
void FixNVEBodyAgent::determine_next_reneighbor()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (length(avec->bonus[atom->body[i]].dvalue) > L_max) {
        next_reneighbor = update->ntimestep;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  maxtag_all = current max atom ID for all atoms
------------------------------------------------------------------------- */

void FixNVEBodyAgent::find_maxid()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
  find how many added atoms in this timestep, across all processors
------------------------------------------------------------------------- */

// void FixNVEBodyAgent::find_nadded_atoms()
// {
//   tagint *tag = atom->tag;
//   tagint *molecule = atom->molecule;
//   int nlocal = atom->nlocal;

//   tagint max = 0;
//   for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
//   MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
// }