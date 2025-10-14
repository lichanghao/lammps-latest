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