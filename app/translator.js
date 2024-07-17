

const { scope } = require("./scope.js");
const { tracker } = require("../parsing/errors.js");
const { TokenType } = require("../parsing/pratt/tokentype.js");
const {
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("../parsing/pratt/expressions.js");
const { FunctionDefinition } = require("./input_expressions.js");

function sortByDependence(lines) {
    return lines.sort((a, b) => {
        const aDependsOnB = a.requirements.includes(b.name);
        const bDependsOnA = b.requirements.includes(a.name);

        if (aDependsOnB && bDependsOnA) {
            tracker.error("circular definitions");
            return 0;
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
    lines = sortByDependence(lines);

    let glslString = "";

    for (const line of lines) {
        let args = "";
        if (line instanceof FunctionDefinition) {
            const locals = (scope.builtin[line.name]) ? Object.keys(scope.builtin[line.name].locals) : Object.keys(scope.userGlobal[line.name].locals);
            args = locals.sort(key => locals[key].index).map(req => "vec2 " + req).join(",");
        }
        const body = astToGLSL(line.ast);
        const alias = scope.builtin[line.name]?.shaderAlias ?? scope.userGlobal[line.name]?.shaderAlias;
        glslString += `vec2 ${alias}(${args}) {\n\treturn ${body};\n}\n\n`;
    }

    return glslString;
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
        const alias = scope.builtin[ast.mFunction]?.shaderAlias ?? scope.userGlobal[ast.mFunction]?.shaderAlias;
        return `${alias}(${ast.mArgs.map(astToGLSL).join(",")})`;
    } else if (ast instanceof NameExpression) {
        const alias = scope.builtin[ast.mName]?.shaderAlias ?? scope.userGlobal[ast.mName]?.shaderAlias;
        // if (!(ast.mName.slice(0, 4) === "udf_") || !newVars.includes(ast.mName)) {
        if (alias === undefined) {
            // it's a local
            return ast.mName;
        } else if (!(alias.slice(0, 4) === "udf_")) {
            return alias;
        } else {
            return `${alias}()`;
        }
    } else if (ast instanceof NumberExpression) {
        return `vec2(${ast.toString()}, 0.)`;
    }
}


module.exports = {
    translateToGLSL,
};