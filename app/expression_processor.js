

/**
 * Goal: classify each expression as a certain type of InputExpression
 * this will make it easier to decide which additional UI components each line needs
 * (such as sliders, rendering options)
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
const {
    FunctionDefinition, VariableDefinition,
    EvaluatableLine,
} = require("./input_expressions.js");