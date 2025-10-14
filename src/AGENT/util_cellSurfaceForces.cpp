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
#include "util_cellSurfaceForces.h"

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "util_integrals.h"

// Expose math functions for both standard and boost MPFR types
using std::sqrt;
using std::abs;

/**
 * Compute the force and moment vectors due to cell-surface repulsion on 
 * the given cell.
 *
 * This function assumes that the z-orientation (n[2]) is nonnegative.  
 *
 * @param r Cell center. 
 * @param n Cell orientation.
 * @param half_l Cell half-length.  
 * @param R Cell radius.
 * @param E0 Elastic modulus of EPS.
 * @returns A six-dimensional vector that contains the force and moment 
 *          vectors due to cell-surface repulsion. 
 */
std::vector<double> cellSurfaceRepulsionForce(std::vector<double>& r,
                                              std::vector<double>& n, 
                                              const double half_l,
                                              const double R, const double E0)
{
    std::vector<double> force {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double EPSILON = 1e-8; 
    const double l = 2 * half_l;  

    // If the z-coordinate of the cell's orientation is near zero ... 
    if (abs(n[2]) < EPSILON)
    {
        // In this case, the force is nonzero only if phi > 0
        double phi = R - r[2]; 
        if (phi > 0)
        {
            // Calculate the force, which is only nonzero in the z-direction;
            // the moment is always zero  
            force[2] = 2 * E0 * phi * l;
        }
    }
    // Otherwise ...
    else
    {
        // Compute the force vector 
        double nz2 = n[2] * n[2];
        double ss = (R - r[2]) / n[2];  
        double int1 = integral1<double>(    // Integral of \delta_i(s)
            r[2], n[2], R, half_l, 1.0, ss
        );
        double int2 = integral1<double>(    // Integral of \sqrt{\delta_i(s)}
            r[2], n[2], R, half_l, 0.5, ss
        );
        force[2] = 2 * E0 * ((1.0 - nz2) * int1 + sqrt(R) * nz2 * int2);

        // Compute the moment vector 
        double int3 = integral2<double>(    // Integral of s * \delta_i(s)
            r[2], n[2], R, half_l, 1.0, ss
        );
        double int4 = integral2<double>(    // Integral of s * \sqrt{\delta_i(s)}
            r[2], n[2], R, half_l, 0.5, ss
        );
        std::vector<double> cross {n[1], -n[0], 0.0};  
        force[3] = 2 * E0 * cross[0] * ((1.0 - nz2) * int3 + sqrt(R) * nz2 * int4); 
        force[4] = 2 * E0 * cross[1] * ((1.0 - nz2) * int3 + sqrt(R) * nz2 * int4); 
    }

    return force;
}

/**
 * Compute the force and moment vectors due to cell-surface adhesion on 
 * the given cell. 
 *
 * This function assumes that the z-orientation (n[2]) is nonnegative.  
 *
 * @param r Cell center. 
 * @param n Cell orientation.
 * @param half_l Cell half-length.  
 * @param R Cell radius.
 * @param sigma0 Cell-surface adhesion energy density. 
 * @returns A six-dimensional vector that contains the force and moment 
 *          vectors due to cell-surface adhesion. 
 */
std::vector<double> cellSurfaceAdhesionForce(std::vector<double>& r,
                                             std::vector<double>& n, 
                                             const double half_l,
                                             const double R,
                                             const double sigma0)
{
    std::vector<double> force {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double EPSILON = 1e-8; 
    const double l = 2 * half_l;  

    // If the z-coordinate of the cell's orientation is near zero ... 
    if (abs(n[2]) < EPSILON)
    {
        // In this case, the force is nonzero only if phi > 0
        double phi = R - r[2]; 
        if (phi > 0)
        {
            // Calculate the force, which is only nonzero in the z-direction;
            // the moment is always zero  
            force[2] = -0.5 * sigma0 * sqrt(R) * l / sqrt(phi); 
        }
    }
    // Otherwise ...
    else
    {
        // Compute the force vector 
        double nz2 = n[2] * n[2];
        double ss = (R - r[2]) / n[2];  
        double int1 = integral1<double>(    // Integral of \delta_i^{-1/2}(s)
            r[2], n[2], R, half_l, -0.5, ss
        );
        double term2 = 0;
        if (abs(ss) <= half_l)
            term2 = boost::math::constants::pi<double>() * sigma0 * R * n[2]; 
        force[2] = -0.5 * sigma0 * sqrt(R) * (1.0 - nz2) * int1 - term2; 

        // Compute the moment vector 
        double int2 = integral2<double>(    // Integral of s * \delta_i^{-1/2}(s)
            r[2], n[2], R, half_l, -0.5, ss
        );
        std::vector<double> cross {n[1], -n[0], 0.0};
        double w1 = 0.0, w2 = 0.0;
        if (abs(ss) <= half_l)
        {
            w1 = boost::math::constants::pi<double>() * sigma0 * R * n[2] * ss * cross[0];
            w2 = boost::math::constants::pi<double>() * sigma0 * R * n[2] * ss * cross[1]; 
        }   
        force[3] = -0.5 * sigma0 * sqrt(R) * (1.0 - nz2) * int2 * cross[0] - w1; 
        force[4] = -0.5 * sigma0 * sqrt(R) * (1.0 - nz2) * int2 * cross[1] - w2;  
    }

    return force;
}

/**
 * Compute the force and moment vectors due to cell-surface friction on 
 * the given cell.
 *
 * This function assumes that the z-orientation (n[2]) is nonnegative.  
 *
 * @param r Cell center. 
 * @param n Cell orientation.
 * @param half_l Cell half-length.
 * @param dr Cell velocity.
 * @param omega Cell angular velocity.   
 * @param R Cell radius.
 * @param eta Cell-surface friction coefficient. 
 * @returns A six-dimensional vector that contains the force and moment 
 *          vectors due to cell-surface friction.
 */
std::vector<double> cellSurfaceFrictionForce(std::vector<double>& r,
                                             std::vector<double>& n, 
                                             const double half_l,
                                             std::vector<double>& dr,
                                             std::vector<double>& omega,  
                                             const double R, const double eta)
{
    std::vector<double> force {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double EPSILON = 1e-8;
    std::vector<double> q {dr[0], dr[1], 0.0};
    std::vector<double> m {
        omega[1] * n[2] - omega[2] * n[1], 
        omega[2] * n[0] - omega[0] * n[2],
        0.0
    };
    const double l = 2 * half_l;  

    // If the z-coordinate of the cell's orientation is near zero ... 
    if (abs(n[2]) < EPSILON)
    {
        // In this case, the force is nonzero only if phi > 0
        double phi = R - r[2]; 
        if (phi > 0)
        {
            // Calculate the force ... 
            double area = sqrt(R) * sqrt(phi) * l; 
            force[0] = -eta * q[0] * area / R;
            force[1] = -eta * q[1] * area / R; 

            // ... and the moment vector, which depends on the cross product
            // of the orientation vector (n) with m  
            std::vector<double> cross {0.0, 0.0, n[0] * m[1] - n[1] * m[0]};  
            force[5] = -eta * cross[2] * area * l * l / (12.0 * R);  
        }
    }
    // Otherwise ...
    else
    {
        // Compute the force vector 
        double ss = (R - r[2]) / n[2];  
        double int1 = areaIntegral1<double>(    // Cell-surface contact area 
            r[2], n[2], R, half_l, ss
        );
        double int2 = areaIntegral2<double>(    // Integral of s * a_i(s)
            r[2], n[2], R, half_l, ss
        );
        double int3 = areaIntegral3<double>(    // Integral of s^2 * a_i(s)
            r[2], n[2], R, half_l, ss 
        );
        force[0] = -eta * (q[0] * int1 + m[0] * int2) / R; 
        force[1] = -eta * (q[1] * int1 + m[1] * int2) / R; 

        // Compute the moment vector, which depends on the cross products of 
        // the orientation vector (n) with q and with m 
        std::vector<double> cross1 {-n[2] * q[1], n[2] * q[0], n[0] * q[1] - n[1] * q[0]}; 
        std::vector<double> cross2 {-n[2] * m[1], n[2] * m[0], n[0] * m[1] - n[1] * m[0]};  
        force[3] = -eta * (cross1[0] * int2 + cross2[0] * int3) / R; 
        force[4] = -eta * (cross1[1] * int2 + cross2[1] * int3) / R;
        force[5] = -eta * (cross1[2] * int2 + cross2[2] * int3) / R;  
    }

    return force; 
}