/**
 * Implementations of cell-surface interaction forces.
 *
 * This is an Eigen-free implementation that takes data for a single cell 
 * as input. This is not called by any function within the other header 
 * files, but can be consulted for a more transparent implementation of 
 * the forces. 
 *
 * Authors:
 *     Kee-Myoung Nam
 *
 * Last updated:
 *     4/15/2025
 */

#ifndef BIOFILM_CELL_SURFACE_FORCES_3D_HPP
#define BIOFILM_CELL_SURFACE_FORCES_3D_HPP

#include <vector>

std::vector<double> cellSurfaceRepulsionForce(std::vector<double>& r,
                                              std::vector<double>& n, 
                                              const double half_l,
                                              const double R, const double E0);


std::vector<double> cellSurfaceAdhesionForce(std::vector<double>& r,
                                             std::vector<double>& n, 
                                             const double half_l,
                                             const double R,
                                             const double sigma0);

std::vector<double> cellSurfaceFrictionForce(std::vector<double>& r,
                                             std::vector<double>& n, 
                                             const double half_l,
                                             std::vector<double>& dr,
                                             std::vector<double>& omega,  
                                             const double R, const double eta);


#endif