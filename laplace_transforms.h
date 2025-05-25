#ifndef LAPLACE_TRANSFORMS_H
#define LAPLACE_TRANSFORMS_H

#include <string>
#include <stdexcept> // For exceptions

// Helper for factorial (n!)
double factorial(int n);

namespace Laplace {
    // Function declarations
    std::string transform_constant(double c);
    std::string transform_t_pow_n(int n, double coeff = 1.0);
    std::string transform_exp(double a, double coeff = 1.0);
    std::string transform_sin(double omega, double coeff = 1.0);
    std::string transform_cos(double omega, double coeff = 1.0);
    std::string transform_t_exp(double a, double coeff = 1.0);
    std::string transform_t_sin(double omega, double coeff = 1.0);
    std::string transform_t_cos(double omega, double coeff = 1.0);
    std::string transform_exp_sin(double a, double omega, double coeff = 1.0);
    std::string transform_exp_cos(double a, double omega, double coeff = 1.0);
    std::string transform_sinh(double omega, double coeff = 1.0);
    std::string transform_cosh(double omega, double coeff = 1.0);
    std::string transform_t_sinh(double omega, double coeff = 1.0);
    std::string transform_t_cosh(double omega, double coeff = 1.0);
    std::string transform_exp_sinh(double a, double omega, double coeff = 1.0);
    std::string transform_exp_cosh(double a, double omega, double coeff = 1.0);

    // New function declarations for combined t*exp*trig/hyperbolic forms
    std::string transform_t_exp_sin(double a, double omega, double coeff = 1.0);
    std::string transform_t_exp_cos(double a, double omega, double coeff = 1.0);
    std::string transform_t_exp_sinh(double a, double omega, double coeff = 1.0);
    std::string transform_t_exp_cosh(double a, double omega, double coeff = 1.0);

} // namespace Laplace

#endif // LAPLACE_TRANSFORMS_H
