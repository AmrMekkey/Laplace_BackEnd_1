#include <iostream>
#include <string>
#include <vector>
#include <limits> // Required for std::numeric_limits
#include <stdexcept> // Required for std::runtime_error
#include "C:\Users\WEST\Programming\C lang\Laplace App\include\laplace_transforms.h"
#include "C:\Users\WEST\Programming\C lang\Laplace App\include\parser.h"

int main() {
    std::cout << "Laplace Transform Calculator" << std::endl;

    std::string input_function;
    Parser parser;
    bool continue_solving = true;

    do {
        // Clear previous state for a new calculation
        parser = Parser(); // Re-initialize the parser for a fresh state
        input_function.clear();

        std::cout << "\nEnter a function of 't' (e.g., 3*t^2 + 2*sin(5*t) - 4*e^(-3*t) + 7):\n";
        std::cout << "Supported combined forms: t*exp(a*t)*sin(omega*t), t*exp(a*t)*cos(omega*t), etc.\n";
        std::cout << "Input: ";
        std::getline(std::cin, input_function);

        try {
            std::vector<ParsedTerm> terms = parser.parse_expression(input_function);
            std::string total_laplace_transform = "";

            for (size_t i = 0; i < terms.size(); ++i) {
                ParsedTerm term = terms[i];
                std::string term_laplace_str;

                switch (term.type) {
                    case FunctionType::CONSTANT:
                        term_laplace_str = Laplace::transform_constant(term.coefficient);
                        break;
                    case FunctionType::T_POW_N:
                        if (term.parameters.empty()) throw std::runtime_error("Missing exponent for t");
                        term_laplace_str = Laplace::transform_t_pow_n(static_cast<int>(term.parameters[0]), term.coefficient);
                        break;
                    case FunctionType::SIN:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for sin");
                        term_laplace_str = Laplace::transform_sin(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::COS:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for cos");
                        term_laplace_str = Laplace::transform_cos(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::EXP:
                        if (term.parameters.empty()) throw std::runtime_error("Missing 'a' for exp");
                        term_laplace_str = Laplace::transform_exp(term.parameters[0], term.coefficient);
                        break;
                    // Basic Hyperbolic Functions
                    case FunctionType::SINH:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for sinh");
                        term_laplace_str = Laplace::transform_sinh(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::COSH:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for cosh");
                        term_laplace_str = Laplace::transform_cosh(term.parameters[0], term.coefficient);
                        break;
                    // Compound Function Types (t*func, exp*func)
                    case FunctionType::T_EXP:
                        if (term.parameters.empty()) throw std::runtime_error("Missing 'a' for t*exp");
                        term_laplace_str = Laplace::transform_t_exp(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::T_SIN:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for t*sin");
                        term_laplace_str = Laplace::transform_t_sin(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::T_COS:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for t*cos");
                        term_laplace_str = Laplace::transform_t_cos(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::EXP_SIN:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for exp*sin");
                        term_laplace_str = Laplace::transform_exp_sin(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    case FunctionType::EXP_COS:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for exp*cos");
                        term_laplace_str = Laplace::transform_exp_cos(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    // Compound Hyperbolic Functions
                    case FunctionType::T_SINH:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for t*sinh");
                        term_laplace_str = Laplace::transform_t_sinh(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::T_COSH:
                        if (term.parameters.empty()) throw std::runtime_error("Missing omega for t*cosh");
                        term_laplace_str = Laplace::transform_t_cosh(term.parameters[0], term.coefficient);
                        break;
                    case FunctionType::EXP_SINH:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for exp*sinh");
                        term_laplace_str = Laplace::transform_exp_sinh(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    case FunctionType::EXP_COSH:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for exp*cosh");
                        term_laplace_str = Laplace::transform_exp_cosh(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    // NEW Combined Functions (t*exp*trig/hyperbolic)
                    case FunctionType::T_EXP_SIN:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for t*exp*sin");
                        term_laplace_str = Laplace::transform_t_exp_sin(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    case FunctionType::T_EXP_COS:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for t*exp*cos");
                        term_laplace_str = Laplace::transform_t_exp_cos(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    case FunctionType::T_EXP_SINH:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for t*exp*sinh");
                        term_laplace_str = Laplace::transform_t_exp_sinh(term.parameters[0], term.parameters[1], term.coefficient);
                        break;
                    case FunctionType::T_EXP_COSH:
                        if (term.parameters.size() < 2) throw std::runtime_error("Missing 'a' or 'omega' for t*exp*cosh");
                        term_laplace_str = Laplace::transform_t_exp_cosh(term.parameters[0], term.parameters[1], term.coefficient);
                        break;

                    case FunctionType::UNRECOGNIZED:
                        term_laplace_str = "[ERROR: Unrecognized Term]";
                        break;
                    default:
                        term_laplace_str = "[ERROR: Unknown Function Type]";
                        break;
                }

                if (i == 0) {
                    total_laplace_transform = term_laplace_str;
                } else {
                    // Check if the current term's coefficient is negative to determine if we should add '-' or '+'
                    // Note: The transform functions themselves embed the coefficient's sign.
                    // So, we only need to add '+' if the term is positive.
                    // If the term_laplace_str starts with '-', it implies a negative coefficient.
                    if (!term_laplace_str.empty() && term_laplace_str[0] == '-') {
                        total_laplace_transform += " " + term_laplace_str;
                    } else {
                        total_laplace_transform += " + " + term_laplace_str;
                    }
                }
            }
            std::cout << "L{" << input_function << "} = " << total_laplace_transform << std::endl;

        } catch (const std::runtime_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
        }

        // Ask user if they want to solve another problem
        std::cout << "\nSolve another problem? (y/n): ";
        std::string user_choice;
        std::getline(std::cin, user_choice); // Read the user's choice

        if (user_choice != "y" && user_choice != "Y") {
            continue_solving = false;
        }

    } while (continue_solving);

    std::cout << "Exiting Laplace Transform Calculator. Goodbye!" << std::endl;
    return 0;
}
