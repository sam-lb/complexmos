/**
 * I think per-pixel function evaluation in javascript is fast enough to run in real time; rendering is the bottleneck
 * so JS generates the mesh and coloring can be done with a dynamic property passed as a uniform in reGL
 */



const { valueScope } = require("./scope.js");
const { tracker } = require("../parsing/errors.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("../parsing/pratt/expressions.js");
const { TokenType } = require("../parsing/pratt/tokentype.js");
const { complex, Complex } = require("../math/complex.js");




function evaluate(ast, args=null) {
    // at this point, the ast has gone through the lexer and parser's checks already
    if (ast instanceof AssignExpression) {
        // the expression is trying to define a function or variable
        const left = ast.mLeft;
        if (left instanceof CallExpression) {
            // the expression is defining a function
            valueScope[left.mFunction.mName] = new Evaluatable(ast.mRight,
                left.mArgs.map((arg) => (arg instanceof OperatorExpression) ? arg.mLeft.mName : arg.mName)
            );
        } else {
            // the expression is defining a variable
            valueScope[left.mName] = new Evaluatable(ast.mRight);
        }
        return null;
    } else {
        return new Evaluatable(ast, args);
    }
}


class Evaluatable {

    constructor(ast, args=null) {
        this.ast = ast;
        this.args = (args === null) ? [] : args;
    }

    call(args) {
        for (const arg of this.args) {
            if (args[arg] === undefined) {
                tracker.error(`function requires argument ${arg}`);
                return null;
            }
        }
        try {
            return this._call(this.ast, args);
        } catch (error) {
            tracker.error(`An unrecognized error has occurred: ${error.message}`);
            return null;
        }
    }

    _call(ast, args=null) {
        args = (args === null) ? {} : args;
        if (ast instanceof OperatorExpression) {
            const arg1 = this._call(ast.mLeft, args);
            const arg2 = this._call(ast.mRight, args);

            switch(ast.mOperator) {
                case (TokenType.PLUS):
                    return Complex.add(arg1, arg2);        
                case (TokenType.MINUS):
                    return Complex.sub(arg1, arg2);
                case (TokenType.ASTERISK):
                    return Complex.mult(arg1, arg2);
                case (TokenType.SLASH):
                    return Complex.div(arg1, arg2);
                case (TokenType.CARET):
                    return Complex.pow(arg1, arg2);
                default:
                    tracker.error("unexpected operator encountered");
                    return null;
            }
        } else if (ast instanceof PrefixExpression) {
            const arg1 = this._call(ast.mRight, args);
            if (ast.mOperator === TokenType.MINUS) {
                return arg1.scale(-1);
            } else {
                tracker.error("unexpected prefix operator encountered");
            }
        } else if (ast instanceof CallExpression) {
            const functionArgs = ast.mArgs.map(arg => this._call(arg, args));
            // we don't check local argument scope here because higher order functions aren't allowed
            if (valueScope[ast.mFunction] !== undefined) {
                const func = valueScope[ast.mFunction];
                if (func instanceof Evaluatable) {
                    const argMap = {};
                    func.args.forEach((key, index) => argMap[key] = functionArgs[index]);
                    return func.call(argMap);
                } else {
                    return func(...functionArgs);
                }
            } else {
                tracker.error(`could not resolve function ${ast.mFunction}. Note: higher order functions are not yet supported`);
            }
        } else if (ast instanceof NameExpression) {
            // console.log(ast.mName, "rear", args);
            if (args[ast.mName] !== undefined) {
                return args[ast.mName];
            } else if (valueScope[ast.mName] !== undefined) {
                let value = valueScope[ast.mName];
                if (value instanceof Evaluatable) value = value.call();
                return value;
            } else {
                tracker.error(`could not resolve variable ${ast.mName}`);
                return null;
            }
        } else if (ast instanceof NumberExpression) {
            return complex(ast.mNumber, 0);
        } else {
            tracker.error("what is even going on");
        }
    }
}



module.exports = {
    evaluate, Evaluatable,
};