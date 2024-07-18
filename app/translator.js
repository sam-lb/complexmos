

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
    // "depends on" does not define a transitive ordering, so this can't be done with a custom compareFn to Array.prototype.sort
    const result = [lines.pop()];
    const dependsOn = (a, b) => a.fullDependence.includes(b.name); // a depends on b

    for (const line of lines) {
        let added = false;
        for (let i=result.length-1; i>=0; i--) {
            const current = result[i];
            if (dependsOn(line, current)) {
                result.splice(i+1, 0, line);
                added = true;
                break;
            }
        }
        if (!added) {
            result.unshift(line);
        }
    }

    return result;
}

function buildFullDependence(lines) {
    const getByName = (name) => lines.filter(line => line.name === name)[0];

    for (const line of lines) {
        let seen = [];
        let needToSearch = [line.name];
        while (needToSearch.length > 0) {
            const dep = getByName(needToSearch.pop());
            needToSearch = needToSearch.concat(dep.requirements.filter(req => !seen.includes(req)));
            seen = seen.concat(dep.requirements);
        }
        line.fullDependence = seen.slice();
    }
}

function checkNoCycles(lines) {
    const dependsOn = (a, b) => a.fullDependence.includes(b.name); // a depends on b

    for (let i=0; i<lines.length; i++) {
        for (let j=i+1; j<lines.length; j++) {
            const line1 = lines[i], line2 = lines[j];
            if (dependsOn(line1, line2) && dependsOn(line2, line1)) {
                tracker.error(`circular defintions for ${line1.name} and ${line2.name}`);
                return;
            }
        }
    }
}

function translateToGLSL(lines) {
    if (lines.length === 0) return;
    buildFullDependence(lines);
    if (tracker.hasError) return;
    checkNoCycles(lines);
    if (tracker.hasError) return;

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