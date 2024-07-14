

const { Lexer } = require("../parsing/pratt/lexer.js");
const { ExpressionParser } = require("../parsing/pratt/expression_parser.js");
const { scope } = require("./scope.js");
const { tracker } = require("../parsing/errors.js");
const { cleanLatex } = require("../parsing/latex_convert.js");
const { TokenType } = require("../parsing/pratt/tokentype.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("../parsing/pratt/expressions.js");
const { FunctionDefinition, VariableDefinition } = require("./input_expressions.js");

let newVars = [];

function sortByDependence(lines) {
    return lines.sort((a, b) => {
        const aDependsOnB = a.requirements.includes(b.name);
        const bDependsOnA = b.requirements.includes(a.name);

        if (aDependsOnB && bDependsOnA) {
            tracker.error("circular definitions");
        } else if (aDependsOnB) {
            return 1;
        } else if (bDependsOnA) {
            return -1;
        } else {
            return 0;
        }
    });
}

function translateToGLSL(lines) {
    if (!(lines instanceof Array)) return {};
    lines = sortByDependence(lines);

    let glslString = "";

    for (const line of lines) {
        let args = "";
        if (line instanceof FunctionDefinition) {
            const locals = (scope.builtin[line.name]) ? Object.keys(scope.builtin[line.name].locals) : Object.keys(scope.userGlobal[line.name].locals);
            args = locals.sort(key => locals[key].index).map(req => "vec2 " + req).join(",");
        }
        const body = astToGLSL(line.ast);
        glslString += `vec2 ${line.name}(${args}) {\n\treturn ${body};\n}\n\n`;
    }

    return glslString;
}

function htranslateToGLSL(fields) {
    newVars = [];
    const expressions = [];
    scope.userGlobal = {"f": {isFunction: true}};
    
    tracker.clear();
    tracker.setCallback((message, target) => console.log(message));
    for (const id of Object.keys(fields)) {
        const latex = cleanLatex(fields[id].field.latex());
        if (!latex.includes("=")) continue; // only assignments are relevant
        expressions.push(latex);
    }

    let tokenArrays = [];
    const lexer = new Lexer(null, true);
    lexer.setScope(scope);
    for (const expr of expressions) {
        lexer.setText(expr);
        lexer.tokenize();
        const tokens = lexer.getTokens();
        const originalName = tokens[0].text;
        for (const token of tokens) {
            if (token.mtype === TokenType.NAME) {
                if (scope.builtin[token.text]?.shaderAlias) {
                    token.text = scope.builtin[token.text].shaderAlias;
                } else {
                    token.text = "udf_" + token.text; // to avoid naming conflicts with glsl's builtins
                }
            }
        }
        tokenArrays.push({
            "name": tokens[0].text,
            "tokens": tokens,
            "dependencies": ((L) => L.slice(1, L.length))(tokens.filter((token) => token.mtype === TokenType.NAME)).map(token => token.text),
        });
        newVars.push(tokens[0].text);
        scope.userGlobal[originalName] = {isFunction: tokens[1].text === "("};
    }

    lexer.setAllowUnboundIdentifiers(false);
    lexer.setScope(scope);
    for (const expr of expressions) {
        lexer.setText(expr);
        lexer.tokenize(); // to accumulate errors
    }
    if ((tracker.hasError) && tracker.message.includes("Undefined")) {
        return { "glsl": "", "valid": false };
    }
    tracker.clear();

    tokenArrays = tokenArrays.sort((a, b) => {
        const aDependsOnB = a.dependencies.includes(b.name);
        const bDependsOnA = b.dependencies.includes(a.name);

        if (aDependsOnB && bDependsOnA) {
            tracker.error("circular definitions");
        } else if (aDependsOnB) {
            return 1;
        } else if (bDependsOnA) {
            return -1;
        } else {
            return 0;
        }
    });

    let glslString = "";

    for (const tokenArray of tokenArrays) {
        const parser = new ExpressionParser(tokenArray.tokens);
        tokenArray.ast = parser.parseExpression();
        if (tracker.hasError) {
            return { "glsl": "", "valid": false };
        }
        
        if (tokenArray.ast.mLeft instanceof CallExpression) {
            const args = [];
            for (const arg of tokenArray.ast.mLeft.mArgs) {
                if (arg instanceof OperatorExpression) {
                    args.push("vec2 " + arg.mLeft.mName);
                } else {
                    args.push("vec2 " + arg.mName);
                }
            }
            const functionBody = astToGLSL(tokenArray.ast.mRight);
            glslString += `vec2 ${tokenArray.ast.mLeft.mFunction}(${args.join(",")}) {\n\treturn ${functionBody};\n}\n\n`;
        } else {
            const definitionBody = astToGLSL(tokenArray.ast.mRight);
            glslString += `vec2 ${tokenArray.ast.mLeft.mName}() {\n\treturn ${definitionBody};\n}\n\n`;
            // constants are translated to functions too because global variables have to be constant expressions in glsl
        }
    }

    if (tokenArrays.some((tokenArr => {
        const left = tokenArr.ast.mLeft;
        if (!(left instanceof CallExpression)) return false;
        if (left.mFunction.mName !== "udf_f") return false;
        if (left.mArgs.length !== 1) return false;
        
        const arg = left.mArgs[0];
        if (arg instanceof OperatorExpression) {
            return arg.mLeft.mName === "z";
        } else {
            return arg.mName === "z";
        }
    })) && !tracker.hasError) {
        return {
            "glsl": glslString,
            "valid": true,
        };
    } else {
        return {
            "glsl": "",
            "valid": false,
        };
    }
}

function astToGLSL(ast) {
    if (ast instanceof OperatorExpression) {
        const arg1 = ast.mLeft;
        const arg2 = ast.mRight;

        switch(ast.mOperator) {
            case (TokenType.PLUS):
                return `addC(${astToGLSL(arg1)}, ${astToGLSL(arg2)})`;
            case (TokenType.MINUS):
                return `subC(${astToGLSL(arg1)}, ${astToGLSL(arg2)})`;
            case (TokenType.ASTERISK):
                return `multC(${astToGLSL(arg1)}, ${astToGLSL(arg2)})`;
            case (TokenType.SLASH):
                return `divC(${astToGLSL(arg1)}, ${astToGLSL(arg2)})`;
            case (TokenType.CARET):
                return `powC(${astToGLSL(arg1)}, ${astToGLSL(arg2)})`;
        }
    } else if (ast instanceof PrefixExpression) {
        return `-${astToGLSL(ast.mRight)}`;
    } else if (ast instanceof CallExpression) {
        return `${ast.mFunction}(${ast.mArgs.map(astToGLSL).join(",")})`;
    } else if (ast instanceof NameExpression) {
        if (!(ast.mName.slice(0, 4) === "udf_") || !newVars.includes(ast.mName)) {
            return ast.mName;
        } else {
            return `${ast.mName}()`;
        }
    } else if (ast instanceof NumberExpression) {
        return `vec2(${ast.toString()}, 0.)`;
    }
}


module.exports = {
    translateToGLSL,
};