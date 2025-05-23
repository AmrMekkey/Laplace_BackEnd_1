#include "C:\Users\WEST\Programming\C lang\Laplace App\include\parser.h"
#include <cctype>               // For isdigit, isalpha, isspace, isalnum
#include <cmath>                // For M_PI (if used directly, but parser.h provides PI_CONSTANT)
#include <sstream>              // For std::stod conversion and general string manipulation
#include <iostream>             // For potential debugging output (can be removed later)


// --- Token Constructors Definition ---
Token::Token(TokenType t, std::string txt) : type(t), text(std::move(txt)), value(0.0) {
    if (type == TokenType::NUMBER && !text.empty()) {
        try {
            value = std::stod(text);
        } catch (const std::out_of_range& oor) {
            throw std::runtime_error("Number out of range: " + text);
        } catch (const std::invalid_argument& ia) {
            // This can happen if text is not a valid double, though tokenizer should ensure it is.
            throw std::runtime_error("Invalid number format: " + text);
        }
    }
}

// Constructor for tokens that are not numbers or have pre-defined text
Token::Token(TokenType t, const char* txt_char) : type(t), text(txt_char), value(0.0) {}

// --- Tokenizer Function Definition ---
std::vector<Token> tokenize(const std::string& input) {
    std::vector<Token> tokens;
    size_t pos = 0;

    while (pos < input.length()) {
        char current_char = input[pos];

        if (std::isspace(current_char)) {
            pos++;
            continue;
        }

        // Handle numbers: allows for leading '.', e.g., ".5"
        if (std::isdigit(current_char) || (current_char == '.' && pos + 1 < input.length() && std::isdigit(input[pos+1]))) {
            std::string num_str;
            bool decimal_found = false;
            while (pos < input.length() && (std::isdigit(input[pos]) || input[pos] == '.')) {
                if (input[pos] == '.') {
                    if (decimal_found) {
                        // Allow only one decimal point
                        throw std::runtime_error("Multiple decimal points in number: " + num_str + input[pos]);
                    }
                    decimal_found = true;
                }
                num_str += input[pos];
                pos++;
            }
            tokens.emplace_back(TokenType::NUMBER, num_str);
            continue;
        }

        // Handle identifiers (function names, 't', 'PI', 'e', 'exp')
        if (std::isalpha(current_char)) {
            std::string ident_str;
            ident_str += current_char;
            pos++;
            while (pos < input.length() && (std::isalnum(input[pos]))) { // allows letters and numbers
                ident_str += input[pos];
                pos++;
            }
            // All identifiers are parsed as IDENTIFIER type, their specific meaning (PI, sin, t)
            // will be interpreted by the parser based on their text.
            tokens.emplace_back(TokenType::IDENTIFIER, ident_str);
            continue;
        }

        // Handle operators and parentheses
        switch (current_char) {
            case '+': tokens.emplace_back(TokenType::PLUS, "+"); pos++; break;
            case '-': tokens.emplace_back(TokenType::MINUS, "-"); pos++; break;
            case '*': tokens.emplace_back(TokenType::MULTIPLY, "*"); pos++; break;
            case '/': tokens.emplace_back(TokenType::DIVIDE, "/"); pos++; break;
            case '^': tokens.emplace_back(TokenType::POWER, "^"); pos++; break;
            case '(': tokens.emplace_back(TokenType::LPAREN, "("); pos++; break;
            case ')': tokens.emplace_back(TokenType::RPAREN, ")"); pos++; break;
            default:
                throw std::runtime_error("Unknown character in input: " + std::string(1, current_char));
                // If you prefer to skip unknown characters rather than throwing:
                // tokens.emplace_back(TokenType::UNKNOWN, std::string(1, current_char));
                // pos++;
                // break;
        }
    }
    tokens.emplace_back(TokenType::END_OF_INPUT, ""); // Mark the end of input
    return tokens;
}

// --- Parser Class Method Definitions ---

Token& Parser::current_token() {
    if (token_idx_ >= tokens_.size()) {
        throw std::runtime_error("Unexpected end of input.");
    }
    return tokens_[token_idx_];
}

Token& Parser::peek_token(size_t offset) {
    if (token_idx_ + offset >= tokens_.size()) {
        throw std::runtime_error("Unexpected end of input on peek.");
    }
    return tokens_[token_idx_ + offset];
}

void Parser::consume_token() {
    if (token_idx_ < tokens_.size()) {
        token_idx_++;
    }
}

// Helper to evaluate simple parameter expressions like "2", "PI", "-2", "2*PI"
// This is primarily for the *coefficient* of a standalone constant,
// or for very simple arguments like 'a' in exp(a*t) or 'omega' in sin(omega*t)
// It assumes parameters are constants or PI. 't' is handled explicitly in argument parsing.
double Parser::evaluate_simple_parameter() {
    double val = 1.0;
    double sign = 1.0;

    if (current_token().type == TokenType::MINUS) {
        sign = -1.0;
        consume_token();
    } else if (current_token().type == TokenType::PLUS) {
        consume_token();
    }

    if (current_token().type == TokenType::NUMBER) {
        val = current_token().value;
        consume_token();
    } else if (current_token().text == "PI") {
        val = PI_CONSTANT;
        consume_token();
    } else {
        throw std::runtime_error("Expected number or PI for parameter value, got: " + current_token().text);
    }

    val *= sign;

    // Handle simple multiplication like "2*PI" or "PI*2" (though usually written 2*PI)
    if (current_token().type == TokenType::MULTIPLY) {
        consume_token(); // Consume '*'
        if (current_token().text == "PI") {
            val *= PI_CONSTANT;
            consume_token();
        } else if (current_token().type == TokenType::NUMBER) {
            val *= current_token().value;
            consume_token();
        } else {
            throw std::runtime_error("Expected PI or number after '*' in parameter, got: " + current_token().text);
        }
    }
    return val;
}

// Helper to parse arguments like "a*t", "omega*t", "a", "omega"
// Returns the constant factor 'a' or 'omega'. Assumes 't' is the variable.
double Parser::evaluate_simple_parameter_argument() {
    double val = 1.0; // Default coefficient if only 't' is present, e.g., sin(t) -> omega=1
    double sign = 1.0;
    bool t_seen = false;

    // Handle leading sign (e.g., exp(-3*t))
    if (current_token().type == TokenType::MINUS) {
        sign = -1.0;
        consume_token();
    } else if (current_token().type == TokenType::PLUS) {
        consume_token(); // consume optional '+'
    }

    // Expected patterns:
    // 1. [NUMBER | PI] [*] [t]
    // 2. [t] [*] [NUMBER | PI] (less common but possible)
    // 3. Just [NUMBER | PI]
    // 4. Just [t]

    // Parse the first component (number, PI, or t)
    if (current_token().type == TokenType::NUMBER) {
        val = current_token().value;
        consume_token();
    } else if (current_token().text == "PI") {
        val = PI_CONSTANT;
        consume_token();
    } else if (current_token().text == "t") {
        val = 1.0; // Implies 1*t
        t_seen = true;
        consume_token();
    } else {
        throw std::runtime_error("Invalid function argument. Expected number, PI, or 't', got: " + current_token().text);
    }

    // Check for multiplication (e.g., 2*t, PI*t, t*2, t*PI)
    if (current_token().type == TokenType::MULTIPLY) {
        consume_token(); // Consume '*'
        if (current_token().type == TokenType::NUMBER) {
            if (t_seen) { // If 't*NUMBER', update coefficient
                val *= current_token().value;
            } else { // If 'NUMBER*NUMBER' or 'PI*NUMBER', update coefficient
                val *= current_token().value;
            }
            consume_token();
        } else if (current_token().text == "PI") {
            if (t_seen) { // If 't*PI', update coefficient
                val *= PI_CONSTANT;
            } else { // If 'NUMBER*PI' or 'PI*PI', update coefficient
                val *= PI_CONSTANT;
            }
            consume_token();
        } else if (current_token().text == "t") {
            if (t_seen) { // Already saw 't', like 't*t' (unsupported)
                throw std::runtime_error("Unexpected 't' in function argument (e.g., 't*t' is not supported).");
            }
            t_seen = true;
            consume_token();
        } else {
            throw std::runtime_error("Expected number, PI, or 't' after '*' in argument, got: " + current_token().text);
        }
    }

    // Final check for 't' if not seen yet (e.g., '2t' - which is NOT supported by your assumption
    // but this check ensures 't' is present if the function expects a 't' dependent parameter).
    // Given your "multiplication is explicit" assumption, `2t` won't be parsed.
    // This `t_seen` check is mostly for ensuring we correctly interpret 't' in `sin(t)` vs `sin(2)`.
    // For Laplace functions, we assume arguments are `a*t` or `omega*t`.
    // If a function like sin(2) is given, it's considered sin(2*t) with omega = 2.
    // If it's a constant like sin(2) (not sin(2*t)), it should be treated as a CONSTANT term,
    // but the current parser's design for function arguments will try to extract 'omega' or 'a'.
    // We are proceeding with the assumption that arguments for sin/cos/exp always relate to 't'
    // in the form of `omega*t` or `a*t`, where 'omega'/'a' can be 1 if just 't' is present.

    return val * sign;
}


// Tries to parse a single term which might have a coefficient
// Helper function (private member of Parser class) to parse a single simple function (t, sin, cos, exp)
// and return a temporary ParsedTerm for it. This makes the main parse_term_with_coeff cleaner.
// This function should be placed inside the Parser class definition in .h or just before its usage in .cpp
// We'll define it as a private helper function of Parser.
ParsedTerm Parser::parse_single_function_factor() {
    ParsedTerm temp_factor;
    temp_factor.coefficient = 1.0; // Factors initially have coeff 1, main coeff applied later.

    if (current_token().text == "t") {
        consume_token(); // Consume "t"
        temp_factor.type = FunctionType::T_POW_N;
        if (current_token().type == TokenType::POWER) {
            consume_token(); // Consume '^'
            if (current_token().type == TokenType::NUMBER) {
                temp_factor.parameters.push_back(current_token().value); // n
                consume_token();
            } else {
                throw std::runtime_error("Expected exponent (number) after '^' for t^n.");
            }
        } else {
            temp_factor.parameters.push_back(1.0); // Default to t^1 if no power specified
        }
    } else if (current_token().text == "sin") {
        consume_token(); // Consume "sin"
        temp_factor.type = FunctionType::SIN;
        if (current_token().type != TokenType::LPAREN) throw std::runtime_error("Expected '(' after sin.");
        consume_token(); // Consume '('
        temp_factor.parameters.push_back(evaluate_simple_parameter_argument()); // Parse omega
        if (current_token().type != TokenType::RPAREN) throw std::runtime_error("Expected ')' after sin argument.");
        consume_token(); // Consume ')'
    } else if (current_token().text == "cos") {
        consume_token(); // Consume "cos"
        temp_factor.type = FunctionType::COS;
        if (current_token().type != TokenType::LPAREN) throw std::runtime_error("Expected '(' after cos.");
        consume_token(); // Consume '('
        temp_factor.parameters.push_back(evaluate_simple_parameter_argument()); // Parse omega
        if (current_token().type != TokenType::RPAREN) throw std::runtime_error("Expected ')' after cos argument.");
        consume_token(); // Consume ')'
    } else if (current_token().text == "exp" || current_token().text == "e") {
        std::string func_name = current_token().text;
        consume_token(); // Consume "exp" or "e"
        temp_factor.type = FunctionType::EXP;
        if (func_name == "e" && current_token().type == TokenType::POWER) { // For "e^(expr)"
            consume_token(); // Consume '^'
        }
        if (current_token().type != TokenType::LPAREN) throw std::runtime_error("Expected '(' after " + func_name + (func_name == "e" ? "^" : "") + ".");
        consume_token(); // Consume '('
        temp_factor.parameters.push_back(evaluate_simple_parameter_argument()); // Parse 'a'
        if (current_token().type != TokenType::RPAREN) throw std::runtime_error("Expected ')' after exp/e argument.");
        consume_token(); // Consume ')'
    }else if (current_token().text == "sinh") { // NEW: Handle sinh
        consume_token(); // Consume "sinh"
        temp_factor.type = FunctionType::SINH;
        if (current_token().type != TokenType::LPAREN) throw std::runtime_error("Expected '(' after sinh.");
        consume_token(); // Consume '('
        temp_factor.parameters.push_back(evaluate_simple_parameter_argument()); // Parse omega
        if (current_token().type != TokenType::RPAREN) throw std::runtime_error("Expected ')' after sinh argument.");
        consume_token(); // Consume ')'
    } else if (current_token().text == "cosh") { // NEW: Handle cosh
        consume_token(); // Consume "cosh"
        temp_factor.type = FunctionType::COSH;
        if (current_token().type != TokenType::LPAREN) throw std::runtime_error("Expected '(' after cosh.");
        consume_token(); // Consume '('
        temp_factor.parameters.push_back(evaluate_simple_parameter_argument()); // Parse omega
        if (current_token().type != TokenType::RPAREN) throw std::runtime_error("Expected ')' after cosh argument.");
        consume_token(); // Consume ')'
    } else {
        throw std::runtime_error("Expected a function (t, sin, cos, exp, sinh, cosh) but got: " + current_token().text);
    }
    return temp_factor;
}


// REMOVE THE OLD parse_term_with_coeff AND REPLACE IT WITH THIS NEW ONE
void Parser::parse_term_with_coeff(double overall_sign) {
    ParsedTerm final_term;
    double current_coeff = overall_sign; // Start with the sign determined by '+' or '-'

    // Handle explicit leading coefficient (e.g., "3*sin(t)", "5")
    if (current_token().type == TokenType::NUMBER) {
        current_coeff *= current_token().value;
        consume_token();
        if (current_token().type != TokenType::MULTIPLY) {
            // It's a standalone constant term (e.g., "5", "-3")
            final_term.type = FunctionType::CONSTANT;
            final_term.coefficient = current_coeff;
            parsed_terms_.push_back(final_term);
            return;
        }
        consume_token(); // Consume '*'
    }
    final_term.coefficient = current_coeff; // Apply the parsed coefficient to the final term

    // Parse the first function/factor
    ParsedTerm first_factor = parse_single_function_factor();
    final_term.type = first_factor.type;
    final_term.parameters = first_factor.parameters;

    // Check for multiplication with another function/factor
    if (current_token().type == TokenType::MULTIPLY) {
        consume_token(); // Consume '*'

        // Parse the second function/factor
        ParsedTerm second_factor = parse_single_function_factor();

        // Now, combine first_factor and second_factor into a single compound function
        // Check for (t * exp), (t * sin), (t * cos)
        // Check for (t * exp), (t * sin), (t * cos), (t * sinh), (t * cosh)
        if (first_factor.type == FunctionType::T_POW_N && first_factor.parameters[0] == 1.0) { // Check if it's 't' (t^1)
            if (second_factor.type == FunctionType::EXP) {
                final_term.type = FunctionType::T_EXP;
                final_term.parameters = second_factor.parameters; // 'a' from exp(at)
            } else if (second_factor.type == FunctionType::SIN) {
                final_term.type = FunctionType::T_SIN;
                final_term.parameters = second_factor.parameters; // 'omega' from sin(omega*t)
            } else if (second_factor.type == FunctionType::COS) {
                final_term.type = FunctionType::T_COS;
                final_term.parameters = second_factor.parameters; // 'omega' from cos(omega*t)
            } else if (second_factor.type == FunctionType::SINH) { // NEW: t * sinh
                final_term.type = FunctionType::T_SINH;
                final_term.parameters = second_factor.parameters; // 'omega' from sinh(omega*t)
            } else if (second_factor.type == FunctionType::COSH) { // NEW: t * cosh
                final_term.type = FunctionType::T_COSH;
                final_term.parameters = second_factor.parameters; // 'omega' from cosh(omega*t)
            } else {
                throw std::runtime_error("Unsupported multiplication of 't' with " + second_factor.text_representation());
            }
        }
        // Check for (exp * sin), (exp * cos), (exp * sinh), (exp * cosh)
        else if (first_factor.type == FunctionType::EXP) {
            if (second_factor.type == FunctionType::SIN) {
                final_term.type = FunctionType::EXP_SIN;
                final_term.parameters.push_back(first_factor.parameters[0]);  // 'a' from exp(at)
                final_term.parameters.push_back(second_factor.parameters[0]); // 'omega' from sin(omega*t)
            } else if (second_factor.type == FunctionType::COS) {
                final_term.type = FunctionType::EXP_COS;
                final_term.parameters.push_back(first_factor.parameters[0]);  // 'a' from exp(at)
                final_term.parameters.push_back(second_factor.parameters[0]); // 'omega' from cos(omega*t)
            } else if (second_factor.type == FunctionType::SINH) { // NEW: exp * sinh
                final_term.type = FunctionType::EXP_SINH;
                final_term.parameters.push_back(first_factor.parameters[0]);  // 'a' from exp(at)
                final_term.parameters.push_back(second_factor.parameters[0]); // 'omega' from sinh(omega*t)
            } else if (second_factor.type == FunctionType::COSH) { // NEW: exp * cosh
                final_term.type = FunctionType::EXP_COSH;
                final_term.parameters.push_back(first_factor.parameters[0]);  // 'a' from exp(at)
                final_term.parameters.push_back(second_factor.parameters[0]); // 'omega' from cosh(omega*t)
            } else {
                throw std::runtime_error("Unsupported multiplication of exp(at) with " + second_factor.text_representation());
            }
        }
        
        // Handle cases like "t^2 * cos(...)"
        else if (first_factor.type == FunctionType::T_POW_N && first_factor.parameters[0] > 1.0) {
            // For now, we'll keep this as an error for simplicity.
            // Supporting L{t^n * f(t)} requires a different Laplace transform rule or a more advanced parser.
            // If you implement L{t^n * f(t)} = (-1)^n * d^n/ds^n F(s), you'd need the parser to just
            // give you f(t) and n, then main.cpp would call L{f(t)} and differentiate.
            throw std::runtime_error("Multiplication with t^n (n>1) is not directly supported: " + first_factor.text_representation() + " * " + second_factor.text_representation());
        }
        else {
            throw std::runtime_error("Unsupported multiplication of functions: " + first_factor.text_representation() + " * " + second_factor.text_representation());
        }
    }

    parsed_terms_.push_back(final_term); // Add the parsed term (simple or compound)
}

// Add this new helper function's declaration to Parser class in parser.h (as a private member)
// Add this definition (ParsedTerm Parser::parse_single_function_factor()) to parser.cpp
// just before the updated Parser::parse_term_with_coeff.


std::vector<ParsedTerm> Parser::parse_expression(const std::string& input) {
    tokens_ = tokenize(input); // Tokenize the input string
    token_idx_ = 0;             // Reset token index for new parse
    parsed_terms_.clear();      // Clear any previous parsed terms

    if (tokens_.empty() || tokens_[0].type == TokenType::END_OF_INPUT) {
        return parsed_terms_; // Empty input, return empty list
    }

    // Handle leading sign for the first term (e.g., "-t^2 + 5")
    double overall_sign = 1.0;
    if (current_token().type == TokenType::PLUS) {
        consume_token();
    } else if (current_token().type == TokenType::MINUS) {
        overall_sign = -1.0;
        consume_token();
    }

    parse_term_with_coeff(overall_sign); // Parse the first term

    // Loop to parse subsequent terms separated by '+' or '-'
    while (current_token().type == TokenType::PLUS || current_token().type == TokenType::MINUS) {
        overall_sign = (current_token().type == TokenType::PLUS) ? 1.0 : -1.0;
        consume_token(); // Consume PLUS or MINUS operator
        parse_term_with_coeff(overall_sign); // Parse the next term
    }

    // After parsing all terms, ensure we are at the end of input
    if (current_token().type != TokenType::END_OF_INPUT) {
        throw std::runtime_error("Unexpected token at end of expression: " + current_token().text);
    }
    return parsed_terms_;
}