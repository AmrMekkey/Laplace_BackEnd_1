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
        // For Gamma function (which generalizes factorial), tgamma(n+1) can be used.
        // However, standard Laplace transform tables use integer n for t^n.
        throw std::invalid_argument("Factorial is not defined for negative integers.");
    }
    if (n > 20) { // Factorials grow very quickly, watch for overflow with double
        // Consider using tgamma(n + 1) from <cmath> for larger n if precision allows
        // or indicate that result might be too large or symbolic 'n!' is better
    }
    double res = 1.0;
    for (int i = 1; i <= n; ++i) {
        res *= i;
    }
    return res;
}

// Namespace for Laplace Transform functions
namespace Laplace {

/**
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

    if (total_coeff == 0.0 && fact_n != 0.0) { // e.g. coeff is extremely small
         oss << "0";
         return oss.str();
    }
    if (total_coeff != 1.0 || n == 0) { // if n=0, fact_n=1, so print coeff unless it's 1/s
         oss << total_coeff;
    } else if (fact_n == 0 && n > 0) { // Should not happen with positive n and factorial logic
         oss << "0"; // Or handle error
         return oss.str();
    }


    if (n == 0) { // L{c*t^0} = L{c} = c/s
        oss << "/s";
    } else {
        if (total_coeff != 1.0 && fact_n != 0) oss << "/s^" << (n + 1);
        else if (fact_n != 0) oss << "1/s^" << (n+1); // handles t^n case where coeff is 1
        else if (total_coeff != 0) oss << "/s^" << (n+1); // handles coeff * t^n when fact_n might be 1
        else oss.str("0"); // All factors led to zero
    }
    // A refinement to avoid "1/s^power" and just show "1/s^power"
    std::string temp_str = oss.str();
    if (temp_str.rfind("1/", 0) == 0 && total_coeff == 1.0 && n > 0) { // starts with "1/"
        // This logic is a bit tricky; ensure it correctly simplifies.
        // e.g. if total_coeff is 1, it might output "1/s^2". If it's just "t", it should be "1/s^2".
        // if it is "1*t", it is "1/s^2".
        // The main goal is to avoid 1*X -> X, unless X starts with number.
    } else if (total_coeff == 1.0 && n == 0) { // For L{1}, it should be "1/s"
        return "1/s";
    } else if (total_coeff == 0.0) {
        return "0";
    }


    // Simplification for output "1.0*X" to "X" or "Y/s^N"
    if (oss.str().empty() && total_coeff != 0.0) { // case t^0
        oss << total_coeff << "/s";
    } else if (total_coeff == 1.0 && n > 0) {
        oss.str(""); // Clear
        oss << fact_n << "/s^" << (n+1); // This might be redundant if fact_n is 1
        if(fact_n == 1.0) oss.str("1/s^" + std::to_string(n+1));
    } else if (total_coeff != 1.0 && total_coeff != 0.0) {
        // Handled by initial oss << total_coeff
    }


    std::string result_str = oss.str();
    // Final check for "1*X" or "1.0*X"
    if (result_str.rfind(std::to_string(total_coeff),0) == 0 && total_coeff == 1.0 && n > 0) {
         return result_str.substr(std::to_string(total_coeff).length());
    }
    if (result_str.empty() && coeff != 0) return transform_constant(coeff); // for t^0
    if (result_str.empty() && coeff == 0) return "0";

    // Re-evaluating the string construction for clarity
    oss.str(""); oss.clear(); // Clear stream state
    if (total_coeff == 1.0 && n > 0) { // e.g. t, t^2
        oss << "1/s^" << (n + 1);
    } else if (n == 0) { // constant term from t^0
        oss << total_coeff << "/s";
    }
    else {
        oss << total_coeff << "/s^" << (n + 1);
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
    oss << -coeff;
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
    oss << -coeff * 2.0 * omega << "*s";
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
    if (coeff != 1.0) oss << coeff << "*";
    oss << "("<<omega * omega << " - s^2" <<")";
    oss << "/( (s^2 + " << omega * omega << ")^2 )";
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
         if (a != 0.0) {
            if (a > 0.0) oss << " - " << a;
            else oss << " + " << -a;
        }
        oss << ")";
    } else { // coeff is 1
        oss << "(s";
        if (a != 0.0) {
            if (a > 0.0) oss << " - " << a;
            else oss << " + " << -a;
        }
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

// ... (your existing transform_exp_cos function) ...

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
    oss << coeff * -2.0 * omega << "*s";
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
    if (coeff != 1.0) oss << -coeff << "*";
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


}
