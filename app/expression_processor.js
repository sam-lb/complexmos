

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
const { ExpressionParser } = require("../parsing/pratt/expression_parser.js");
const { complex, Complex } = require("../math/complex.js");
const {
    FunctionDefinition, VariableDefinition,
    EvaluatableLine,
} = require("./input_expressions.js");
const { cleanLatex } = require("../parsing/latex_convert.js");


function populateGlobalUserScope(fields) {
    const lexer = new Lexer(null, true);
    // intentionally NOT setting the lexer's scope here
    const functionAssignments = {};

    for (const id in fields) {
        const field = fields[id];
        const latex = cleanLatex(field.field.latex());

        if (!latex.includes("=")) continue;

        lexer.setText(latex);
        lexer.tokenize();
        if (tracker.hasError) return null;
        const tokens = lexer.getTokens();
        const name = tokens[0];

        if (!(name.mtype === TokenType.NAME)) {
            tracker.error("Invalid assignment: left hand side must begin with identifier");
            return null;
        }

        if (!(scope.builtin[name.text] === undefined)) {
            tracker.error(`Cannot overwrite builtin ${name.text}`);
            return null;
        }

        const second = tokens[1]; // not undefined since there's at least an identifier and = at this point
        if (second.mtype === TokenType.ASSIGN) {
            scope.userGlobal[name.text] = {
                isFunction: false,
            };
        } else if (second.mtype === TokenType.ASTERISK && tokens[2]?.mtype === TokenType.LEFT_PAREN) {
            // account for implicit multiplication
            tokens.splice(1, 1);
            scope.userGlobal[name.text] = {
                isFunction: true,
            };
            functionAssignments[name.text] = tokens;
        } else {
            tracker.error("Invalid assignment: left hand side must be identifier or function with argument list");
            return null;
        }
    }

    return functionAssignments;
}

function populateLocalUserScopes(functionAssignments) {
    // specify local variables for functions
    for (const name in functionAssignments) {
        const assignment = [];
        for (const token of functionAssignments[name]) {
            if (token.text === "=") break;
            assignment.push(token);
        }
        const locals = {};
        console.log(assignment, "assignment");
        const ast = (new ExpressionParser(assignment)).parseExpression();
        if (tracker.hasError) return null;
        // if (!(ast instanceof AssignExpression) || !(ast.mLeft instanceof CallExpression)) {
        if (!(ast instanceof CallExpression)) {
            tracker.error("Invalid assignment"); // not a lot of detail in the error message because it's not clear when this might happen
            return null;
        }

        const args = ast.mArgs;
        for (const arg of args) {
            let argName;
            if (arg instanceof NameExpression) {
                // argument without type spec
                argName = arg.mName;
            } else if (arg instanceof OperatorExpression && arg.mOperator === TokenType.COLON) {
                // argument with type spec
                // ignore type for now, since only valid type is complex
                argName = arg.mLeft.mName;
            } else {
                tracker.error(`Invalid argument ${arg.toString()}`);
                return null;
            }

            if (Object.keys(scope.builtin).includes(argName) && !scope.builtin[argName].isParameter) {
                tracker.error(`Cannot locally overwrite builtin identifier ${argName}`);
                return null;
            } else if (Object.keys(scope.userGlobal).includes(argName)) {
                tracker.error(`Cannot locally overwrite globally defined identifier ${argName}`);
                return null;
            }

            locals[argName] = {
                isFunction: false,
                type: "complex",
            };
        }

        scope.userGlobal[name]["locals"] = locals;
    }
    return 1;
}


function populateUserScope(fields) {
    scope.userGlobal = {};
    tracker.clear();

    tracker.setCallback((message, target) => console.log(message));
    const functionAssignments = populateGlobalUserScope(fields);
    if (functionAssignments === null) return;
    const success = populateLocalUserScopes(functionAssignments);
    if (success === null) return;
}


function classifyInput(fields) {
    tracker.setCallback((message, target) => console.log(message));
    const lexer = new Lexer(null, false);
    lexer.setScope(scope);
    const inputExpressions = {
        "functions": [],
        "variables": [],
        "evaluatables": [],
    };

    for (const id in fields) {
        const field = fields[id];
        const latex = cleanLatex(field.field.latex());
        if (latex === "") {
            // skip empty lines
            continue;
        }
        lexer.setText(latex);
        lexer.tokenize();
        const tokens = lexer.getTokens();
        
        if (latex.includes("=")) {
            if (tokens[1]?.mtype === TokenType.LEFT_PAREN) {
                // re-tokenize with local scope
                lexer.setText(latex);
                lexer.setLocalScope(scope.userGlobal[tokens[0].text].locals);
                lexer.tokenize();
                inputExpressions["functions"].push(new FunctionDefinition(lexer.getTokens(), id));
            } else {
                inputExpressions["variables"].push(new VariableDefinition(tokens, id));
            }
        } else {
            inputExpressions["evaluatables"].push(new EvaluatableLine(tokens, id));
        }
    }

    return inputExpressions;
}

function allRequirementsSatisfied(lines, names) {
    if (lines.some(line => {
        for (const req of line.requirements) {
            if (!names.includes(req)) {
                tracker.error(`Unbound variable ${req}`);
                return true;
            }
        }
        return false;
    })) {
        return false;
    }

    return true;
}

function noInvalidRequirements(varsAndFuncs, lines) {
    // check that there are no circular requirements (including self requirements)
    // check that there are no repeated definitions
    for (const line of varsAndFuncs) {
        for (const req of line.requirements) {
            const definitions = lines.filter(def => def.name === req);
            if (definitions.length > 1) {
                tracker.error(`Multiple defintions found for ${req}`);
                return false;
            }
            const definition = definitions[0];
            if (definition.requirements.includes(line.name)) {
                tracker.error(`Circular definition: ${req}, ${line.name}`);
                return false;
            }
        }
    }

    return true;
}

function noRepeatDefinitions(names) {
    console.log(names);
    if (!(Array.from(new Set(names)).length === names.length)) {
        tracker.error("Repeated definitions");
        return false;
    }
    return true;
}


function validateLines(lines) {
    const varsAndFuncs = lines["functions"].concat(lines["variables"]);
    const allLines = Array.prototype.concat(lines["functions"], lines["variables"], lines["evaluatables"]);
    const names = varsAndFuncs.map(line => line.name);

    if (!allRequirementsSatisfied(allLines, names)) return false;
    if (!noInvalidRequirements(varsAndFuncs, allLines)) return false;
    if (!noRepeatDefinitions(names)) return false;

    return true;
}


function validateAST(ast) {
    if (
        ast instanceof AssignExpression ||
        ast instanceof OperatorExpression
    ) {
        validateAST(ast.mLeft);
        if (tracker.hasError) return;
        validateAST(ast.mRight);
    } else if (ast instanceof PrefixExpression) {
        validateAST(ast.mRight);
    } else if (ast instanceof CallExpression) {
        // Validate argument count (and later, type);
        let requiredArgCount;
        if (scope.userGlobal[ast.mFunction]?.locals) {
            requiredArgCount = Object.keys(scope.userGlobal[ast.mFunction].locals).length;
        } else {
            requiredArgCount = Object.keys(scope.builtin[ast.mFunction].locals).length;
        }
        const passedArgCount = ast.mArgs.length;
        if (requiredArgCount !== passedArgCount) {
            tracker.error(`Wrong number of arguments passed to ${ast.mFunction} (expected ${requiredArgCount}, received ${passedArgCount})`);
            return;
        }
        for (const arg of ast.mArgs) {
            validateAST(arg);
            if (tracker.hasError) return;
        }
    }
}


module.exports = {
    classifyInput,
    validateLines,
    populateUserScope,
    validateAST,
};