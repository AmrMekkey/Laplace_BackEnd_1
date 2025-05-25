#include "C:\Users\WEST\Programming\C lang\Laplace App\include\laplace_transforms.h"
#include <string>
#include <vector>
#include <stdexcept> // For exceptions
#include <cmath>     // For std::pow, std::tgamma (for factorial approximation if needed)
#include <sstream>   // For string formatting
#include <iomanip>   // For std::fixed, std::setprecision
#include <iostream>

// Define PI if not available or for consistency
const double PI_CONST = 3.14159265358979323846;

// Helper for factorial (n!)
double factorial(int n) {
    if (n < 0) {
        // Factorial is not defined for negative numbers in this context.
        throw std::invalid_argument("Factorial is not defined for negative integers.");
    }
    if (n > 20) { // Factorials grow very quickly, watch for overflow with double
        // For larger n, consider using tgamma(n + 1) from <cmath> if precision allows,
        // or indicate that the result might be too large or symbolic 'n!' is better.
        // For this application, n is typically small for t^n.
    }
    double res = 1.0;
    for (int i = 1; i <= n; ++i) {
        res *= i;
    }
    return res;
}

// Namespace for Laplace Transform functions
namespace Laplace {

/*
 * @brief Computes the Laplace Transform of a constant c.
 * L{c} = c/s
 * @param c The constant value.
 * @return String representation of the transform.
 */
std::string transform_constant(double c) {
    if (c == 0.0) return "0";
    std::ostringstream oss;
    oss << c << "/s";
    return oss.str();
}

/**
 * @brief Computes the Laplace Transform of t^n.
 * L{coeff * t^n} = coeff * n! / s^(n+1)
 * @param n The exponent (non-negative integer).
 * @param coeff The coefficient multiplying t^n (default is 1.0).
 * @return String representation of the transform.
 */
std::string transform_t_pow_n(int n, double coeff) {
    if (coeff == 0.0) return "0";
    if (n < 0) throw std::invalid_argument("n must be a non-negative integer for L{t^n}.");

    std::ostringstream oss;
    double fact_n = factorial(n);
    double total_coeff = coeff * fact_n;

    if (total_coeff == 0.0) {
        return "0";
    }

    // Simplification logic for 1/s^n+1 vs coeff/s^n+1
    if (n == 0) { // L{coeff*t^0} = L{coeff} = coeff/s
        oss << total_coeff << "/s";
    } else {
        if (total_coeff == 1.0) {
            oss << "1/s^" << (n + 1);
        } else {
            oss << total_coeff << "/s^" << (n + 1);
        }
    }
    return oss.str();
}


/**
 * @brief Computes the Laplace Transform of e^(at).
 * L{coeff * e^(at)} = coeff / (s - a)
 * @param a The constant in the exponent.
 * @param coeff The coefficient multiplying e^(at) (default is 1.0).
 * @return String representation of the transform.
 */
std::string transform_exp(double a, double coeff) {
    if (coeff == 0.0) return "0";
    std::ostringstream oss;
    oss << coeff;
    oss << "/(s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")";
    return oss.str();
}

/**
 * @brief Computes the Laplace Transform of sin(omega*t).
 * L{coeff * sin(omega*t)} = coeff * omega / (s^2 + omega^2)
 * @param omega The angular frequency.
 * @param coeff The coefficient multiplying sin(omega*t) (default is 1.0).
 * @return String representation of the transform.
 */
std::string transform_sin(double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // sin(0) = 0

    std::ostringstream oss;
    oss << coeff * omega;
    oss << "/(s^2 + " << omega * omega << ")";
    return oss.str();
}

/**
 * @brief Computes the Laplace Transform of cos(omega*t).
 * L{coeff * cos(omega*t)} = coeff * s / (s^2 + omega^2)
 * @param omega The angular frequency.
 * @param coeff The coefficient multiplying cos(omega*t) (default is 1.0).
 * @return String representation of the transform.
 */
std::string transform_cos(double omega, double coeff ) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return transform_constant(coeff); // cos(0) = 1, so L{coeff*1}

    std::ostringstream oss;
    if (coeff != 1.0) oss << coeff << "*";
    oss << "s";
    oss << "/(s^2 + " << omega * omega << ")";
    return oss.str();
}

// --- More Advanced/Combined Forms (using properties) ---

/**
 * @brief L{coeff * t * e^(at)} = coeff / (s-a)^2
 */
std::string transform_t_exp(double a, double coeff ) {
    if (coeff == 0.0) return "0";
    std::ostringstream oss;
    oss << coeff; // The formula is 1/(s-a)^2, so coeff goes to numerator
    oss << "/((s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * sin(omega*t)} = coeff * 2*omega*s / (s^2 + omega^2)^2
 */
std::string transform_t_sin(double omega, double coeff ) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // t*sin(0) = 0
    std::ostringstream oss;
    // L{t*sin(wt)} = -d/ds(w/(s^2+w^2)) = - (0 - w*2s) / (s^2+w^2)^2 = 2ws / (s^2+w^2)^2
    oss << coeff * 2.0 * omega << "*s";
    oss << "/((s^2 + " << omega * omega << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * cos(omega*t)} = coeff * (s^2 - omega^2) / (s^2 + omega^2)^2
 */
std::string transform_t_cos(double omega, double coeff ) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return transform_t_pow_n(1, coeff); // t*cos(0) = t

    std::ostringstream oss;
    // L{t*cos(wt)} = -d/ds(s/(s^2+w^2)) = - (1*(s^2+w^2) - s*2s) / (s^2+w^2)^2
    // = - (s^2+w^2 - 2s^2) / (s^2+w^2)^2 = - (w^2 - s^2) / (s^2+w^2)^2 = (s^2 - w^2) / (s^2+w^2)^2
    if (coeff != 1.0) oss << coeff << "*";
    oss << "(s^2 - " << omega * omega << ")";
    oss << "/((s^2 + " << omega * omega << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * e^(at) * sin(omega*t)} = coeff * omega / ((s-a)^2 + omega^2)
 */
std::string transform_exp_sin(double a, double omega, double coeff ) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // e^(at)*sin(0) = 0

    std::ostringstream oss;
    oss << coeff * omega;
    oss << "/((s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")^2 + " << (omega * omega) << ")";
    return oss.str();
}

/**
 * @brief L{coeff * e^(at) * cos(omega*t)} = coeff * (s-a) / ((s-a)^2 + omega^2)
 */
std::string transform_exp_cos(double a, double omega, double coeff ) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return transform_exp(a, coeff); // e^(at)*cos(0) = e^(at)

    std::ostringstream oss;
    if (coeff != 1.0) { // If coeff is 1, it will be (s-a)/...
        oss << coeff << "*(s";
        if (a > 0.0) oss << " - " << a;
        else if (a < 0.0) oss << " + " << -a;
        oss << ")";
    } else { // coeff is 1
        oss << "(s";
        if (a > 0.0) oss << " - " << a;
        else if (a < 0.0) oss << " + " << -a;
        oss << ")";
    }

    oss << "/((s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")^2 + " << (omega * omega) << ")";
    return oss.str();
}

/**
 * @brief Computes the Laplace Transform of sinh(omega*t).
 * L{coeff * sinh(omega*t)} = coeff * omega / (s^2 - omega^2)
 * @param omega The angular frequency.
 * @param coeff The coefficient multiplying sinh(omega*t).
 * @return String representation of the transform.
 */
std::string transform_sinh(double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // sinh(0) = 0

    std::ostringstream oss;
    oss << coeff * omega;
    oss << "/(s^2 - " << omega * omega << ")";
    return oss.str();
}

/**
 * @brief Computes the Laplace Transform of cosh(omega*t).
 * L{coeff * cosh(omega*t)} = coeff * s / (s^2 - omega^2)
 * @param omega The angular frequency.
 * @param coeff The coefficient multiplying cosh(omega*t).
 * @return String representation of the transform.
 */
std::string transform_cosh(double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return Laplace::transform_constant(coeff); // cosh(0) = 1, so L{coeff*1}

    std::ostringstream oss;
    if (coeff != 1.0) oss << coeff << "*";
    oss << "s";
    oss << "/(s^2 - " << omega * omega << ")";
    return oss.str();
}

/**
 * @brief L{coeff * t * sinh(omega*t)} = coeff * 2*omega*s / (s^2 - omega^2)^2
 * @param omega The angular frequency.
 * @param coeff The coefficient.
 * @return String representation of the transform.
 */
std::string transform_t_sinh(double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // t*sinh(0) = 0

    std::ostringstream oss;
    // L{t*sinh(wt)} = -d/ds(w/(s^2-w^2)) = - (0 - w*2s) / (s^2-w^2)^2 = 2ws / (s^2-w^2)^2
    oss << coeff * 2.0 * omega << "*s";
    oss << "/((s^2 - " << omega * omega << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * cosh(omega*t)} = coeff * (s^2 + omega^2) / (s^2 - omega^2)^2
 * @param omega The angular frequency.
 * @param coeff The coefficient.
 * @return String representation of the transform.
 */
std::string transform_t_cosh(double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return Laplace::transform_t_pow_n(1, coeff); // t*cosh(0) = t

    std::ostringstream oss;
    // L{t*cosh(wt)} = -d/ds(s/(s^2-w^2)) = - (1*(s^2-w^2) - s*2s) / (s^2-w^2)^2
    // = - (s^2-w^2 - 2s^2) / (s^2-w^2)^2 = - (-w^2 - s^2) / (s^2-w^2)^2 = (s^2 + w^2) / (s^2-w^2)^2
    if (coeff != 1.0) oss << coeff << "*";
    oss << "(s^2 + " << omega * omega << ")";
    oss << "/((s^2 - " << omega * omega << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * e^(at) * sinh(omega*t)} = coeff * omega / ((s-a)^2 - omega^2)
 * @param a The constant in the exponent.
 * @param omega The angular frequency.
 * @param coeff The coefficient.
 * @return String representation of the transform.
 */
std::string transform_exp_sinh(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // e^(at)*sinh(0) = 0

    std::ostringstream oss;
    oss << coeff * omega;
    oss << "/((s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")^2 - " << (omega * omega) << ")";
    return oss.str();
}

/**
 * @brief L{coeff * e^(at) * cosh(omega*t)} = coeff * (s-a) / ((s-a)^2 - omega^2)
 * @param a The constant in the exponent.
 * @param omega The angular frequency.
 * @param coeff The coefficient.
 * @return String representation of the transform.
 */
std::string transform_exp_cosh(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return Laplace::transform_exp(a, coeff); // e^(at)*cosh(0) = e^(at)

    std::ostringstream oss;
    if (coeff != 1.0) {
        oss << coeff << "*(s";
        if (a > 0.0) oss << " - " << a;
        else if (a < 0.0) oss << " + " << -a;
        oss << ")";
    } else {
        oss << "(s";
        if (a > 0.0) oss << " - " << a;
        else if (a < 0.0) oss << " + " << -a;
        oss << ")";
    }
    oss << "/((s";
    if (a != 0.0) {
        if (a > 0.0) oss << " - " << a;
        else oss << " + " << -a;
    }
    oss << ")^2 - " << (omega * omega) << ")";
    return oss.str();
}


// --- New Combined Transform Functions ---

/**
 * @brief L{coeff * t * e^(at) * sin(omega*t)} = coeff * 2*omega*(s-a) / (((s-a)^2 + omega^2)^2)
 * Derived from L{t*f(t)} = -d/ds(F(s)) and L{e^(at)f(t)} = F(s-a)
 * F(s) for sin(wt) is w/(s^2+w^2). Shifted F(s-a) is w/((s-a)^2+w^2).
 * -d/ds [w/((s-a)^2+w^2)] = - [0 - w*2(s-a)] / ((s-a)^2+w^2)^2 = 2w(s-a) / ((s-a)^2+w^2)^2
 */
std::string transform_t_exp_sin(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // t*e^(at)*sin(0) = 0

    std::ostringstream oss;
    oss << coeff * 2.0 * omega << "*(s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")";
    oss << "/(((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 + " << (omega * omega) << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * e^(at) * cos(omega*t)} = coeff * ((s-a)^2 - omega^2) / (((s-a)^2 + omega^2)^2)
 * Derived from L{t*f(t)} = -d/ds(F(s)) and L{e^(at)f(t)} = F(s-a)
 * F(s) for cos(wt) is s/(s^2+w^2). Shifted F(s-a) is (s-a)/((s-a)^2+w^2).
 * -d/ds [(s-a)/((s-a)^2+w^2)] = - [1*((s-a)^2+w^2) - (s-a)*2(s-a)] / ((s-a)^2+w^2)^2
 * = - [ (s-a)^2+w^2 - 2(s-a)^2 ] / ((s-a)^2+w^2)^2
 * = - [ w^2 - (s-a)^2 ] / ((s-a)^2+w^2)^2 = ((s-a)^2 - w^2) / ((s-a)^2+w^2)^2
 */
std::string transform_t_exp_cos(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return transform_t_exp(a, coeff); // t*e^(at)*cos(0) = t*e^(at)

    std::ostringstream oss;
    if (coeff != 1.0) oss << coeff << "*";
    oss << "((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 - " << (omega * omega) << ")";
    oss << "/(((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 + " << (omega * omega) << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * e^(at) * sinh(omega*t)} = coeff * 2*omega*(s-a) / (((s-a)^2 - omega^2)^2)
 * Derived from L{t*f(t)} = -d/ds(F(s)) and L{e^(at)f(t)} = F(s-a)
 * F(s) for sinh(wt) is w/(s^2-w^2). Shifted F(s-a) is w/((s-a)^2-w^2).
 * -d/ds [w/((s-a)^2-w^2)] = - [0 - w*2(s-a)] / ((s-a)^2-w^2)^2 = 2w(s-a) / ((s-a)^2-w^2)^2
 */
std::string transform_t_exp_sinh(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return "0"; // t*e^(at)*sinh(0) = 0

    std::ostringstream oss;
    oss << coeff * 2.0 * omega << "*(s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")";
    oss << "/(((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 - " << (omega * omega) << ")^2)";
    return oss.str();
}

/**
 * @brief L{coeff * t * e^(at) * cosh(omega*t)} = coeff * ((s-a)^2 + omega^2) / (((s-a)^2 - omega^2)^2)
 * Derived from L{t*f(t)} = -d/ds(F(s)) and L{e^(at)f(t)} = F(s-a)
 * F(s) for cosh(wt) is s/(s^2-w^2). Shifted F(s-a) is (s-a)/((s-a)^2-w^2).
 * -d/ds [(s-a)/((s-a)^2-w^2)] = - [1*((s-a)^2-w^2) - (s-a)*2(s-a)] / ((s-a)^2-w^2)^2
 * = - [ (s-a)^2-w^2 - 2(s-a)^2 ] / ((s-a)^2-w^2)^2
 * = - [ -w^2 - (s-a)^2 ] / ((s-a)^2-w^2)^2 = ((s-a)^2 + w^2) / ((s-a)^2-w^2)^2
 */
std::string transform_t_exp_cosh(double a, double omega, double coeff) {
    if (coeff == 0.0) return "0";
    if (omega == 0.0) return transform_t_exp(a, coeff); // t*e^(at)*cosh(0) = t*e^(at)

    std::ostringstream oss;
    if (coeff != 1.0) oss << coeff << "*";
    oss << "((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 + " << (omega * omega) << ")";
    oss << "/(((s";
    if (a > 0.0) oss << " - " << a;
    else if (a < 0.0) oss << " + " << -a;
    oss << ")^2 - " << (omega * omega) << ")^2)";
    return oss.str();
}

} // namespace Laplace
