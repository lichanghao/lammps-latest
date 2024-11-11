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

#ifndef LMP_SURF_EXTRA_H
#define LMP_SURF_EXTRA_H

namespace SurfExtra {

  struct Connect2d {      // line connectivity
    int np1,np2;          // # of lines connected to pts 1,2 (including self)
    int *neigh_p1;        // indices of all lines connected to pt1 (if np1 > 1)
    int *neigh_p2;        // ditto for pt2
    int flags;            // future flags for end pt coupling
  };

  struct Connect3d {      // tri connectivity
    int ne1,ne2,ne3;      // # of tris connected to edges 1,2,3 (including self)
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3 (including self)
    int *neigh_e1;        // indices of all tris connected to 1-2 edge (if ne1 > 1)
    int *neigh_e2;        // ditto for 2-3 edge
    int *neigh_e3;        // ditto for 3-1 edge
    int *neigh_c1;        // indices of all tris connected to corner pt 1 (if nc1 > 1)
    int *neigh_c2;        // ditto for corner pt 2
    int *neigh_c3;        // ditto for corner pt 3
    int flags;            // future flags for edge and corner pt coupling
  };

  int overlap_sphere_line(double *, double, double *, double *,
                          double *, double *, double &);
  int overlap_sphere_tri(double *, double,
                         double *, double *, double *, double *,
                         double *, double *, double &);
  int nearest_point_line(double *, double *, double *, double *);

}

#endif
