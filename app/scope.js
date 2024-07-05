
const { complex, Complex } = require("../math/complex.js");


const scope = {
    builtin: {
        "z": {
            isFunction: false,
            shaderAlias: "z",
        },

        "norm": {
            isFunction: true,
            shaderAlias: "normC",
        },
        "normSq": {
            isFunction: true,
            shaderAlias: "normSqC",
        },
        "arg": {
            isFunction: true,
            shaderAlias: "argC",
        },
        "inv": {
            isFunction: true,
            shaderAlias: "invC",
        },
        "exp": { 
            isFunction: true,
            shaderAlias: "expC",
        },
        "ln": { 
            isFunction: true,
            shaderAlias: "lnC",
        },
        "sqrt": { 
            isFunction: true,
            shaderAlias: "sqrtC",
        },
        "sin": { 
            isFunction: true,
            shaderAlias: "sinC",
        },
        "cos": { 
            isFunction: true,
            shaderAlias: "cosC",
        },
        "tan": {
            isFunction: true,
            shaderAlias: "tanC",
        },
        "asin": {
            isFunction: true,
            shaderAlias: "asinC",
        },
        "acos": {
            isFunction: true,
            shaderAlias: "acosC",
        },
        "atan": {
            isFunction: true,
            shaderAlias: "atanC",
        },
        "sinh": {
            isFunction: true,
            shaderAlias: "sinhC",
        },
        "cosh": {
            isFunction: true,
            shaderAlias: "coshC",
        },
        "tanh": {
            isFunction: true,
            shaderAlias: "tanhC",
        },
        "re": {
            isFunction: true,
            shaderAlias: "reC",
        },
        "im": {
            isFunction: true,
            shaderAlias: "imC",
        },
        "Gamma": {
            isFunction: true,
            shaderAlias: "GammaC",
        },
        "beta": {
            isFunction: true,
            shaderAlias: "betaC",
        },
        "min": {
            isFunction: true,
            shaderAlias: "minC",
        },
        "max": {
            isFunction: true,
            shaderAlias: "maxC",
        },
        "lerp": {
            isFunction: true,
            shaderAlias: "lerpC",
        },
        "conj": {
            isFunction: true,
            shaderAlias: "conjC",
        },

        "i": {
            isFunction: false,
            shaderAlias: "i",
        },
        "pi": {
            isFunction: false,
            shaderAlias: "piC",
        },
        "tau": {
            isFunction: false,
            shaderAlias: "tpiC",
        },
        "e": {
            isFunction: false,
            shaderAlias: "e",
        },
        "realBounds": {
            isFunction: false,
            shaderAlias: "xBounds",
        },
        "imagBounds": {
            isFunction: false,
            shaderAlias: "yBounds",
        },

        /* datatypes (for type annotation) */
        "complex": { isFunction: false, isType: true, },
        "array": { isFunction: false, isType: true, },
        "matrix": { isFunction: false, isType: true, },
        "function": { isFunction: false, isType: true, },
    },

    userGlobal: {
    },
}


const defaultValueScope = {
    "z": complex(1, 0),
    "norm": (z) => complex(z.norm(), 0),
    "normSq": (z) => complex(z.normSq(), 0),
    "arg": (z) => complex(z.arg(), 0),
    "inv": Complex.inv,
    "exp": Complex.exp,
    "ln": Complex.ln,
    "sqrt": Complex.sqrt,
    "sin": Complex.sin,
    "cos": Complex.cos,
    "tan": Complex.tan,
    "acos": Complex.acos,
    "asin": Complex.asin,
    "atan": Complex.atan,
    "sinh": Complex.sinh,
    "cosh": Complex.cosh,
    "tanh": Complex.tanh,
    "re": (z) => complex(z.re, 0),
    "im": (z) => complex(z.im, 0),
    "Gamma": Complex.gamma,
    "beta": Complex.beta,
    "min": Complex.min,
    "max": Complex.max,
    "lerp": (z1, z2, t) => Complex.mult(complex(1, 0).sub(t), z1).add(z2.mult(t)),
    "conj": (z) => z.conj(),
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
    "realBounds": complex(0, 0),
    "imagBounds": complex(0, 0),
};

const valueScope = {};

module.exports = {
    scope, defaultValueScope, valueScope,
};