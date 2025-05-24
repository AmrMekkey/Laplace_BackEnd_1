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
        if (std::isdigit(current_char) || ( (current_char == '.')  && pos + 1 < input.length() && std::isdigit(input[pos+1]))) { // Check for leading digit or decimal point
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



// Parses the smallest unit: a number, a variable 't', a function call, or a parenthesized expression.
ParsedTerm Parser::parse_factor() {
    ParsedTerm term;
    term.coefficient = 1.0; // Factors initially have coeff 1.0

    // Handle potential leading sign for a factor itself, e.g., "-(t)" or "-(sin(t))"
    // This is distinct from the overall_sign for addition/subtraction.
    // A more robust parser would handle unary minus here and pass it down.
    // For now, let's assume `evaluate_simple_parameter_argument` handles signs within arguments,
    // and overall term sign is handled by `parse_expression`.
    // If you want to handle "-(t)" as a factor, uncomment and modify this:
    /*
    double factor_sign = 1.0;
    if (current_token().type == TokenType::MINUS) {
        factor_sign = -1.0;
        term.original_term_str += current_token().text;
        consume_token();
    } else if (current_token().type == TokenType::PLUS) {
        term.original_term_str += current_token().text;
        consume_token();
    }
    term.coefficient *= factor_sign;
    */

    // Handle numbers or constants like 'PI'
    if (current_token().type == TokenType::NUMBER) {
        term.type = FunctionType::CONSTANT;
        term.coefficient *= current_token().value; // Apply explicit coefficient
        term.original_term_str += current_token().text;
        consume_token();
        return term;
    } else if (current_token().type == TokenType::IDENTIFIER && current_token().text == "PI") {
        term.type = FunctionType::CONSTANT;
        term.coefficient *= PI_CONSTANT;
        term.original_term_str += current_token().text;
        consume_token();
        return term;
    } else if (current_token().type == TokenType::IDENTIFIER && current_token().text == "t") {
        term.type = FunctionType::T_POW_N;
        term.parameters.push_back(1.0); // Default to t^1
        term.original_term_str += current_token().text;
        consume_token();
        if (current_token().type == TokenType::POWER) { // Handle t^n
            term.original_term_str += current_token().text; // Add '^'
            consume_token(); // Consume '^'
            if (current_token().type != TokenType::NUMBER) {
                throw std::runtime_error("Expected number for exponent after 't^', got: " + current_token().text);
            }
            term.parameters[0] = current_token().value; // Update power 'n'
            term.original_term_str += current_token().text; // Add 'n'
            consume_token();
        }
        return term;
    } else if (current_token().type == TokenType::IDENTIFIER) {
        // It's a function name: sin, cos, exp, sinh, cosh
        std::string func_name = current_token().text;
        term.original_term_str += current_token().text;
        consume_token(); // Consume function name

        if (func_name == "e") { // Special handling for 'e' followed by '^' for exp(at)
            if (current_token().type == TokenType::POWER) {
                term.original_term_str += current_token().text; // Add '^'
                consume_token(); // Consume '^'
            } else {
                // 'e' without '^' or '(', treat as a constant, if it's not a function.
                // For now, assume 'e' implies 'exp' if followed by '(', otherwise error or constant.
                // If you want 'e' as a constant, similar to 'PI', handle it here:
                // term.type = FunctionType::CONSTANT;
                // term.coefficient *= std::exp(1.0); // Mathematical constant e
                // return term;
                throw std::runtime_error("Identifier 'e' must be followed by '^' for exponentiation or '(' for exp() function: " + current_token().text);
            }
        }

        if (current_token().type != TokenType::LPAREN) {
            throw std::runtime_error("Expected '(' after function name " + func_name + ", got: " + current_token().text);
        }
        term.original_term_str += current_token().text; // Add '('
        consume_token(); // Consume '('

        // Parse argument (e.g., 2*t, -3*t, t, PI*t)
        double param_val = evaluate_simple_parameter_argument();
        term.parameters.push_back(param_val);
        // original_term_str for parameter part is built inside evaluate_simple_parameter_argument if needed

        if (current_token().type != TokenType::RPAREN) {
            throw std::runtime_error("Expected ')' after function arguments, got: " + current_token().text);
        }
        term.original_term_str += current_token().text; // Add ')'
        consume_token(); // Consume ')'

        if (func_name == "sin") {
            term.type = FunctionType::SIN;
        } else if (func_name == "cos") {
            term.type = FunctionType::COS;
        } else if (func_name == "exp" || func_name == "e") { // 'e' handled as exp
            term.type = FunctionType::EXP;
        } else if (func_name == "sinh") {
            term.type = FunctionType::SINH;
        } else if (func_name == "cosh") {
            term.type = FunctionType::COSH;
        } else {
            throw std::runtime_error("Unrecognized function name: " + func_name);
        }
        return term;
    } else if (current_token().type == TokenType::LPAREN) {
        term.original_term_str += current_token().text; // Add '('
        consume_token(); // Consume '('
        // Recursively parse the entire sub-expression within parentheses
        // This is a simpler version; if you need full expression parsing inside parens,
        // you'd call parse_expression or a similar top-level rule here.
        // For now, we assume simple terms inside like (5*sin(t))
        std::vector<ParsedTerm> sub_terms = parse_expression_in_parentheses_helper(); // Call helper for recursive parsing
        if (sub_terms.size() != 1) {
            throw std::runtime_error("Parenthesized expressions currently only support a single combined term for Laplace transform purposes.");
        }
        term = sub_terms[0]; // Take the result of the sub-expression as the factor

        if (current_token().type != TokenType::RPAREN) {
            throw std::runtime_error("Expected ')' after parenthesized expression, got: " + current_token().text);
        }
        term.original_term_str += current_token().text; // Add ')'
        consume_token(); // Consume ')'
        return term;
    } else {
        throw std::runtime_error("Unexpected token while parsing factor: " + current_token().text);
    }
}


// Implementation of combine_multiplied_terms
// Make sure this is part of the Parser class (i.e., Parser:: prefix)
ParsedTerm Parser::combine_multiplied_terms(const ParsedTerm& term1, const ParsedTerm& term2) {
    ParsedTerm combined_term;
    combined_term.coefficient = term1.coefficient * term2.coefficient;
    combined_term.original_term_str = "(" + term1.original_term_str + "*" + term2.original_term_str + ")";

    // Case 1: exp(at) * sin(omega t) or sin(omega t) * exp(at)
    if ((term1.type == FunctionType::EXP && term2.type == FunctionType::SIN) ||
        (term1.type == FunctionType::SIN && term2.type == FunctionType::EXP)) {
        if (term1.type == FunctionType::EXP) {
            combined_term.type = FunctionType::EXP_SIN;
            combined_term.parameters.push_back(term1.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term2.parameters[0]); // 'omega' from sin
        } else { // term1 is SIN
            combined_term.type = FunctionType::EXP_SIN;
            combined_term.parameters.push_back(term2.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term1.parameters[0]); // 'omega' from sin
        }
    }
    // Case 2: exp(at) * cos(omega t) or cos(omega t) * exp(at)
    else if ((term1.type == FunctionType::EXP && term2.type == FunctionType::COS) ||
             (term1.type == FunctionType::COS && term2.type == FunctionType::EXP)) {
        if (term1.type == FunctionType::EXP) {
            combined_term.type = FunctionType::EXP_COS;
            combined_term.parameters.push_back(term1.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term2.parameters[0]); // 'omega' from cos
        } else { // term1 is COS
            combined_term.type = FunctionType::EXP_COS;
            combined_term.parameters.push_back(term2.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term1.parameters[0]); // 'omega' from cos
        }
    }
    // Case 3: exp(at) * sinh(omega t) or sinh(omega t) * exp(at)
    else if ((term1.type == FunctionType::EXP && term2.type == FunctionType::SINH) ||
             (term1.type == FunctionType::SINH && term2.type == FunctionType::EXP)) {
        if (term1.type == FunctionType::EXP) {
            combined_term.type = FunctionType::EXP_SINH;
            combined_term.parameters.push_back(term1.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term2.parameters[0]); // 'omega' from sinh
        } else { // term1 is SINH
            combined_term.type = FunctionType::EXP_SINH;
            combined_term.parameters.push_back(term2.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term1.parameters[0]); // 'omega' from sinh
        }
    }
    // Case 4: exp(at) * cosh(omega t) or cosh(omega t) * exp(at)
    else if ((term1.type == FunctionType::EXP && term2.type == FunctionType::COSH) ||
             (term1.type == FunctionType::COSH && term2.type == FunctionType::EXP)) {
        if (term1.type == FunctionType::EXP) {
            combined_term.type = FunctionType::EXP_COSH;
            combined_term.parameters.push_back(term1.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term2.parameters[0]); // 'omega' from cosh
        } else { // term1 is COSH
            combined_term.type = FunctionType::EXP_COSH;
            combined_term.parameters.push_back(term2.parameters[0]); // 'a' from exp
            combined_term.parameters.push_back(term1.parameters[0]); // 'omega' from cosh
        }
    }
    // Case 5: t * sin(omega t) or sin(omega t) * t
    else if ((term1.type == FunctionType::T_POW_N && term1.parameters[0] == 1.0 && term2.type == FunctionType::SIN) ||
             (term1.type == FunctionType::SIN && term2.type == FunctionType::T_POW_N && term2.parameters[0] == 1.0)) {
        combined_term.type = FunctionType::T_SIN;
        combined_term.parameters.push_back((term1.type == FunctionType::SIN) ? term1.parameters[0] : term2.parameters[0]);
    }
    // Case 6: t * cos(omega t) or cos(omega t) * t
    else if ((term1.type == FunctionType::T_POW_N && term1.parameters[0] == 1.0 && term2.type == FunctionType::COS) ||
             (term1.type == FunctionType::COS && term2.type == FunctionType::T_POW_N && term2.parameters[0] == 1.0)) {
        combined_term.type = FunctionType::T_COS;
        combined_term.parameters.push_back((term1.type == FunctionType::COS) ? term1.parameters[0] : term2.parameters[0]);
    }
    // Case 7: t * sinh(omega t) or sinh(omega t) * t
    else if ((term1.type == FunctionType::T_POW_N && term1.parameters[0] == 1.0 && term2.type == FunctionType::SINH) ||
             (term1.type == FunctionType::SINH && term2.type == FunctionType::T_POW_N && term2.parameters[0] == 1.0)) {
        combined_term.type = FunctionType::T_SINH;
        combined_term.parameters.push_back((term1.type == FunctionType::SINH) ? term1.parameters[0] : term2.parameters[0]);
    }
    // Case 8: t * cosh(omega t) or cosh(omega t) * t
    else if ((term1.type == FunctionType::T_POW_N && term1.parameters[0] == 1.0 && term2.type == FunctionType::COSH) ||
             (term1.type == FunctionType::COSH && term2.type == FunctionType::T_POW_N && term2.parameters[0] == 1.0)) {
        combined_term.type = FunctionType::T_COSH;
        combined_term.parameters.push_back((term1.type == FunctionType::COSH) ? term1.parameters[0] : term2.parameters[0]);
    }
    // Case 9: t * exp(at) or exp(at) * t
    else if ((term1.type == FunctionType::T_POW_N && term1.parameters[0] == 1.0 && term2.type == FunctionType::EXP) ||
             (term1.type == FunctionType::EXP && term2.type == FunctionType::T_POW_N && term2.parameters[0] == 1.0)) {
        combined_term.type = FunctionType::T_EXP;
        combined_term.parameters.push_back((term1.type == FunctionType::EXP) ? term1.parameters[0] : term2.parameters[0]);
    }
    // If we multiply a known function by a constant, it's still the same function type
    // e.g., (2) * sin(t) is still SIN, coefficient just gets multiplied
    else if (term1.type == FunctionType::CONSTANT) {
        combined_term.type = term2.type;
        combined_term.parameters = term2.parameters;
    } else if (term2.type == FunctionType::CONSTANT) {
        combined_term.type = term1.type;
        combined_term.parameters = term1.parameters;
    }
    else {
        // Fallback: If not a recognized combined form, treat as UNRECOGNIZED or throw
        combined_term.type = FunctionType::UNRECOGNIZED;
        throw std::runtime_error("Unsupported multiplication of functions: " + term1.original_term_str + " * " + term2.original_term_str);
    }

    return combined_term;
}

// Parses terms with multiplication and division precedence
ParsedTerm Parser::parse_multiplication() {
    ParsedTerm left_factor = parse_factor(); // Get the first factor

    while (current_token().type == TokenType::MULTIPLY || current_token().type == TokenType::DIVIDE) {
        TokenType op_type = current_token().type;
        consume_token(); // Consume '*' or '/'

        ParsedTerm right_factor = parse_factor(); // Get the next factor

        if (op_type == TokenType::MULTIPLY) {
            left_factor = combine_multiplied_terms(left_factor, right_factor);
        } else { // TokenType::DIVIDE
            // Handle division if you need to. For Laplace, division is less common directly.
            // You could represent it as multiplication by 1/something if 'something' is a constant.
            // For now, throw an error for unsupported division.
            throw std::runtime_error("Division is not supported for Laplace Transform combinations.");
        }
    }
    return left_factor;
}


// Helper function to parse an expression that is expected to be within parentheses.
// It's like a mini-parse_expression, but it expects to find ')' at the end.
std::vector<ParsedTerm> Parser::Parser::parse_expression_in_parentheses_helper() {
    std::vector<ParsedTerm> terms_in_paren;
    double overall_sign_for_term = 1.0;

    // Handle initial sign inside parentheses, e.g. (-sin(t))
    if (current_token().type == TokenType::MINUS) {
        overall_sign_for_term = -1.0;
        consume_token();
    } else if (current_token().type == TokenType::PLUS) {
        consume_token();
    }

    // Parse the first term within parentheses
    ParsedTerm current_parsed_term = parse_multiplication();
    current_parsed_term.coefficient *= overall_sign_for_term;
    terms_in_paren.push_back(current_parsed_term);

    // Continue parsing subsequent terms connected by + or - until RPAREN or END_OF_INPUT
    while (current_token().type == TokenType::PLUS || current_token().type == TokenType::MINUS) {
        overall_sign_for_term = (current_token().type == TokenType::PLUS) ? 1.0 : -1.0;
        consume_token(); // Consume '+' or '-'

        // Parse the next additive term
        current_parsed_term = parse_multiplication();
        current_parsed_term.coefficient *= overall_sign_for_term;
        terms_in_paren.push_back(current_parsed_term);
    }
    // The calling `parse_factor` will consume the RPAREN.
    return terms_in_paren;
}

std::vector<ParsedTerm> Parser::parse_expression(const std::string& input) {
    tokens_ = tokenize(input); // Call the external tokenize function
    token_idx_ = 0; // Reset token index
    // parsed_terms_.clear(); // No longer needed as we return directly

    std::vector<ParsedTerm> result_terms; // Local vector to hold results

    if (tokens_.empty() || current_token().type == TokenType::END_OF_INPUT) {
        return result_terms; // Empty input, return empty list
    }

    // Handle leading sign for the first term
    double overall_sign_for_term = 1.0;
    if (current_token().type == TokenType::MINUS) {
        overall_sign_for_term = -1.0;
        consume_token();
    } else if (current_token().type == TokenType::PLUS) {
        consume_token();
    }

    // Parse the first additive term (which can contain multiplications)
    ParsedTerm current_parsed_term = parse_multiplication(); // Use the new function
    current_parsed_term.coefficient *= overall_sign_for_term;
    result_terms.push_back(current_parsed_term);

    // Continue parsing subsequent additive terms connected by + or -
    while (current_token().type == TokenType::PLUS || current_token().type == TokenType::MINUS) {
        overall_sign_for_term = (current_token().type == TokenType::PLUS) ? 1.0 : -1.0;
        consume_token(); // Consume '+' or '-'

        // Parse the next additive term
        current_parsed_term = parse_multiplication(); // Use the new function
        current_parsed_term.coefficient *= overall_sign_for_term;
        result_terms.push_back(current_parsed_term);
    }

    if (current_token().type != TokenType::END_OF_INPUT) {
        throw std::runtime_error("Unexpected token at end of expression: " + current_token().text);
    }
    return result_terms;
}