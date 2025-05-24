#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>
#include <stdexcept> // For exceptions

// Define PI if not available (M_PI is common in <cmath> but not strictly standard C++)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double PI_CONSTANT = M_PI; // Declare here if needed in header, define in .cpp if only used there.

// --- Tokenizer Types ---
enum class TokenType {
    NUMBER, IDENTIFIER, PLUS, MINUS, MULTIPLY, DIVIDE, POWER, LPAREN, RPAREN, END_OF_INPUT, UNKNOWN,
};

struct Token {
    TokenType type;
    std::string text;
    double value;

    Token(TokenType t, std::string txt = "");
    Token(TokenType t, const char* txt_char);
};

// --- Tokenizer Function Declaration ---
std::vector<Token> tokenize(const std::string& input);

// --- Parser Structures ---
enum class FunctionType {
    UNRECOGNIZED = 0,
    CONSTANT,
    T_POW_N,
    SIN,
    COS,
    EXP,
    SINH,
    COSH,
    T_EXP,
    T_SIN,
    T_COS,
    EXP_SIN,
    EXP_COS,
    T_SINH, 
    T_COSH,    
    EXP_SINH,  
    EXP_COSH,
    UNKNOWN_COMPOUND   
};

struct ParsedTerm {
    double coefficient = 1.0; // Includes sign
    FunctionType type = FunctionType::UNRECOGNIZED;
    std::vector<double> parameters; // For T_POW_N: {n}, EXP: {a}, SIN/COS: {omega}
    std::string original_term_str;

    // // For debugging or easier interpretation
    // std::string original_term_string; // Optional, useful for debugging

    std::string text_representation() const {
        switch (type) {
            case FunctionType::CONSTANT: return "constant";
            case FunctionType::T_POW_N: return "t^" + std::to_string(static_cast<int>(parameters[0]));
            case FunctionType::SIN: return "sin(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::COS: return "cos(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::EXP: return "exp(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::SINH: return "sinh(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::COSH: return "cosh(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::T_EXP: return "t*exp(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::T_SIN: return "t*sin(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::T_COS: return "t*cos(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::EXP_SIN: return "exp(" + std::to_string(parameters[0]) + "*t)*sin(" + std::to_string(parameters[1]) + "*t)";
            case FunctionType::EXP_COS: return "exp(" + std::to_string(parameters[0]) + "*t)*cos(" + std::to_string(parameters[1]) + "*t)";
            case FunctionType::T_SINH: return "t*sinh(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::T_COSH: return "t*cosh(" + std::to_string(parameters[0]) + "*t)";
            case FunctionType::EXP_SINH: return "exp(" + std::to_string(parameters[0]) + "*t)*sinh(" + std::to_string(parameters[1]) + "*t)";
            case FunctionType::EXP_COSH: return "exp(" + std::to_string(parameters[0]) + "*t)*cosh(" + std::to_string(parameters[1]) + "*t)";
            default: return "unrecognized_function";
        }
    }
};

// --- Parser Class Declaration ---
class Parser {
public:
    std::vector<ParsedTerm> parse_expression(const std::string& input);

private:
    std::vector<Token> tokens_;
    size_t token_idx_;
    std::vector<ParsedTerm> parsed_terms_;
    std::vector<ParsedTerm> parse_expression_in_parentheses_helper();

    Token& current_token();
    Token& peek_token(size_t offset = 1);
    void consume_token();
    // double evaluate_simple_parameter();
    // void parse_term_with_coeff(double overall_sign);
    double evaluate_simple_parameter_argument();
    ParsedTerm parse_factor();
    ParsedTerm parse_multiplication(); // Handles terms connected by * or / (higher precedence)
    ParsedTerm combine_multiplied_terms(const ParsedTerm& term1, const ParsedTerm& term2);
};

#endif // PARSER_H