

/**
 * Goal: classify each expression as a certain type of InputExpression
 * this will make it easier to decide which additional UI components each line needs
 * (such as sliders, rendering options)
 */


const { valueScope, scope } = require("./scope.js");
const { tracker } = require("../parsing/errors.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("../parsing/pratt/expressions.js");
const { TokenType } = require("../parsing/pratt/tokentype.js");
const { Lexer } = require("../parsing/pratt/lexer.js");
const { complex, Complex } = require("../math/complex.js");
const {
    FunctionDefinition, VariableDefinition,
    EvaluatableLine,
} = require("./input_expressions.js");
const { cleanLatex } = require("../parsing/latex_convert.js");



function classifyInput(fields) {
    tracker.setCallback((message, target) => console.log(message));
    const lexer = new Lexer(null, true);
    lexer.setScope(scope);
    const inputExpressions = {
        "functions": [],
        "variables": [],
        "evaluatables": [],
    };

    for (const id in fields) {
        const field = fields[id];
        const latex = cleanLatex(field.field.latex());
        lexer.setText(latex);
        lexer.tokenize();
        const tokens = lexer.getTokens();
        
        if (latex.includes("=")) {
            if (tokens[1]?.mtype === TokenType.LEFT_PAREN) {
                inputExpressions["functions"].push(new FunctionDefinition(tokens, id));
            } else {
                inputExpressions["variables"].push(new VariableDefinition(tokens, id));
            }
        } else {
            inputExpressions["evaluatables"].push(new EvaluatableLine(tokens, id));
        }
    }

    console.log(inputExpressions);
    return inputExpressions;
}

function validateLines(lines) {
    // check that none of the function's locals are defined elsewhere
    // check that there are no circular requirements (including self requirements)

    const varsAndFuncs = lines["functions"].concat(lines["variables"]);
    const names = varsAndFuncs.map(line => line.name);

    // check that all the requirements are satisfied
    if (varsAndFuncs.some(line => {
        for (const req of line.requirements) {
            if (!names.includes(req)) {
                tracker.error(`Unbound variable ${req}`);
                return true;
            }
            return false;
        }
    })) {
        return false;
    }

    return true;
}


module.exports = {
    classifyInput,
    validateLines,
};