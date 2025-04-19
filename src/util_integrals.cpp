/**
 * Functions for computing a series of auxiliary integrals involved in 
 * the 3-D biofilm simulations. 
 *
 * Authors:
 *     Kee-Myoung Nam
 *
 * Last updated:
 *     4/16/2025
 */

#ifndef BIOFILM_3D_AUXILIARY_INTEGRALS_HPP
#define BIOFILM_3D_AUXILIARY_INTEGRALS_HPP

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/mpfr.hpp>

// Expose math functions for both standard and boost MPFR types
using std::pow;
using boost::multiprecision::pow;

/**
 * Return the cell-body coordinate at which the cell-surface overlap is zero.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be nonzero.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @returns Cell-body coordinate at which the cell-surface overlap is zero.
 */
template <typename T>
T sstar(const T rz, const T nz, const T R)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_NONZERO
        if (nz == 0)
            throw std::invalid_argument(
                "Cell z-orientation cannot be nonzero when calculating `sstar()`"
            ); 
    #endif
    return (R - rz) / nz;
}

/**
 * Return the cell-surface overlap, \phi_i(s).
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param s Cell-body coordinate. 
 * @returns Cell-surface overlap. 
 */
template <typename T>
T phi(const T rz, const T nz, const T R, const T s)
{
    return R - rz - s * nz;
}

/**
 * Return the cell-surface overlap, \delta_i(s).
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param s Cell-body coordinate. 
 * @returns Cell-surface overlap. 
 */
template <typename T>
T overlap(const T rz, const T nz, const T R, const T s)
{
    T p = phi<T>(rz, nz, R, s);
    if (p > 0)
        return p;
    else
        return 0;
}

/**
 * Return the exponentiated cell-surface overlap, \delta_i^\gamma(s), where 
 * \gamma is a real exponent.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param s Cell-body coordinate. 
 * @param gamma Exponent.
 * @returns Exponentiated cell-surface overlap. 
 */
template <typename T>
T overlapGamma(const T rz, const T nz, const T R, const T s, const T gamma)
{
    T p = phi<T>(rz, nz, R, s); 
    if (p > 0)
        return pow(p, gamma);
    else
        return 0;
}

/**
 * Compute the integral of \delta_i^\gamma(s), where \gamma is a real exponent,
 * from s = -l_i/2 to s = +l_i/2, where l_i is the length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param gamma Exponent; should not be -1.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.  
 */
template <typename T>
T integral1(const T rz, const T nz, const T R, const T half_l, const T gamma,
            const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    if (ss > half_l)
    {
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        T overlap2 = pow(phi<T>(rz, nz, R, half_l), gamma + 1);
        return (overlap1 - overlap2) / (nz * (gamma + 1));
    }
    else if (ss > -half_l)   // -half_l < ss <= half_l
    {
        // overlap2 = 0
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        return overlap1 / (nz * (gamma + 1));
    }
    else                     // ss <= -half_l
    {
        return 0;
    }
}

/**
 * Compute the integral of s * \delta_i^\gamma(s), where \gamma is a real
 * exponent, from s = -l_i/2 to s = +l_i/2, where l_i is the length of the
 * cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius. 
 * @param half_l Half of cell length.
 * @param gamma Exponent; should not be -1 or -2.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.  
 */
template <typename T>
T integral2(const T rz, const T nz, const T R, const T half_l, const T gamma,
            const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    T term1 = integral1<T>(rz, nz, R, half_l, gamma + 1, ss);
    T term2 = 0; 
    if (ss > half_l)
    {
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        T overlap2 = pow(phi<T>(rz, nz, R, half_l), gamma + 1);
        term2 = -half_l * (overlap1 + overlap2); 
    }
    else if (ss > -half_l)    // -half_l < ss <= half_l
    {
        // Overlap at s = ss is zero 
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        term2 = -half_l * overlap1;
    }
    return (term1 + term2) / (nz * (gamma + 1)); 
}

/**
 * Compute the integral of s^2 * \delta_i^\gamma(s), where \gamma is a real
 * exponent, from s = -l_i/2 to s = +l_i/2, where l_i is the length of the
 * cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius. 
 * @param half_l Half of cell length.
 * @param gamma Exponent; should not be -1 or -2.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.  
 */
template <typename T>
T integral3(const T rz, const T nz, const T R, const T half_l, const T gamma,
            const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    T term1 = 2 * integral2<T>(rz, nz, R, half_l, gamma + 1, ss);
    T term2 = 0;
    if (ss > half_l)
    {
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        T overlap2 = pow(phi<T>(rz, nz, R, half_l), gamma + 1);
        term2 = half_l * half_l * (overlap1 - overlap2);
    }
    else if (ss > -half_l)    // -half_l < ss <= half_l
    {
        // Overlap at s = ss is zero 
        T overlap1 = pow(phi<T>(rz, nz, R, -half_l), gamma + 1);
        term2 = half_l * half_l * overlap1;
    }
    return (term1 + term2) / (nz * (gamma + 1));
}

/**
 * Compute the integral of \Theta(\delta_i(s)), where \Theta is the Heaviside
 * step function, from s = -l_i/2 to s = +l_i/2, where l_i is the length of
 * the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`, whose input values are
 *           defined elsewhere and not passed into this function.
 * @returns Desired integral.  
 */
template <typename T>
T integral4(const T half_l, const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    if (ss > half_l)
        return 2 * half_l;    // = half_l - (-half_l)
    else if (ss > -half_l)    // -half_l < ss <= half_l
        return ss + half_l;   // = ss - (-half_l)
    else
        return 0;
}

/**
 * Compute the integral of s * \Theta(\delta_i(s)), where \Theta is the
 * Heaviside step function, from s = -l_i/2 to s = +l_i/2, where l_i is the
 * length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.  
 */
template <typename T>
T integral5(const T rz, const T nz, const T R, const T half_l, const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    if (-half_l < ss && ss <= half_l)
        return (ss * ss - half_l * half_l) / 2.0;
    else 
        return 0.0; 
}

/**
 * Compute the integral of s^2 * \Theta(\delta_i(s)), where \Theta is the
 * Heaviside step function, from s = -l_i/2 to s = +l_i/2, where l_i is the
 * length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.  
 */
template <typename T>
T integral6(const T rz, const T nz, const T R, const T half_l, const T ss)
{
    #ifdef CHECK_CELL_ORIENTATION_ZCOORD_POSITIVE
        if (nz <= 0)
            throw std::invalid_argument(
                "Cell z-orientation must be positive when calculating auxiliary integral"
            ); 
    #endif

    if (ss > half_l)
        return 2 * pow(half_l, 3) / 3.0;
    else if (ss > -half_l)    // -half_l < ss <= half_l
        return (pow(ss, 3) + pow(half_l, 3)) / 3.0; 
    else
        return 0.0;
}

/**
 * Compute the integral of the cell-surface contact area density a_i(s) from
 * s = -l_i/2 to s = +l_i/2, where l_i is the length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.
 */
template <typename T>
T areaIntegral1(const T rz, const T nz, const T R, const T half_l, const T ss)
{
    T term1 = pow(R, 0.5) * (1 - nz * nz) * integral1<T>(rz, nz, R, half_l, 0.5, ss);
    T term2 = boost::math::constants::pi<T>() * R * nz * nz * integral4<T>(half_l, ss);
    return term1 + term2;
}

/**
 * Compute the integral of s * a_i(s), where a_i(s) is the cell-surface
 * contact area density, from s = -l_i/2 to s = +l_i/2, where l_i is the 
 * length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.
 */
template <typename T>
T areaIntegral2(const T rz, const T nz, const T R, const T half_l, const T ss)
{
    T term1 = pow(R, 0.5) * (1 - nz * nz) * integral2<T>(rz, nz, R, half_l, 0.5, ss);
    T term2 = boost::math::constants::pi<T>() * R * nz * nz * integral5<T>(rz, nz, R, half_l, ss);
    return term1 + term2;
}

/**
 * Compute the integral of s^2 * a_i(s), where a_i(s) is the cell-surface
 * contact area density, from s = -l_i/2 to s = +l_i/2, where l_i is the 
 * length of the cell.
 *
 * The z-coordinate of the cell's orientation vector is assumed to be positive.
 *
 * @param rz z-coordinate of cell center.
 * @param nz z-coordinate of cell orientation.
 * @param R Cell radius.
 * @param half_l Half of cell length.
 * @param ss Pre-computed value of `sstar(rz, nz, R)`.
 * @returns Desired integral.
 */
template <typename T>
T areaIntegral3(const T rz, const T nz, const T R, const T half_l, const T ss)
{
    T term1 = pow(R, 0.5) * (1 - nz * nz) * integral3<T>(rz, nz, R, half_l, 0.5, ss);
    T term2 = boost::math::constants::pi<T>() * R * nz * nz * integral6<T>(rz, nz, R, half_l, ss);
    return term1 + term2;
}

#endif
