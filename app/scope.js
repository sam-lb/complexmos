
const { complex, Complex } = require("../math/complex.js");


const scope = {
    builtin: {
        "z": {
            isFunction: false,
            shaderAlias: "z",
            isParameter: true,
        },
        // "t": {
        //     isFunction: false,
        //     shaderAlias: "t",
        //     isParameter: true,
        // },

        "norm": {
            isFunction: true,
            shaderAlias: "normC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "normSq": {
            isFunction: true,
            shaderAlias: "normSqC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "arg": {
            isFunction: true,
            shaderAlias: "argC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "inv": {
            isFunction: true,
            shaderAlias: "invC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "exp": { 
            isFunction: true,
            shaderAlias: "expC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "ln": { 
            isFunction: true,
            shaderAlias: "lnC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "powlog": {
            isFunction: true,
            shaderAlias: "powlogC",
            locals: { "z": { isFunction: false, type: "complex", index: 0, }, "w": { isFunction: false, type: "complex", index: 1, } },
        },
        "sqrt": { 
            isFunction: true,
            shaderAlias: "sqrtC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sin": { 
            isFunction: true,
            shaderAlias: "sinC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "cos": { 
            isFunction: true,
            shaderAlias: "cosC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "tan": {
            isFunction: true,
            shaderAlias: "tanC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "asin": {
            isFunction: true,
            shaderAlias: "asinC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "acos": {
            isFunction: true,
            shaderAlias: "acosC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "atan": {
            isFunction: true,
            shaderAlias: "atanC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sinh": {
            isFunction: true,
            shaderAlias: "sinhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "cosh": {
            isFunction: true,
            shaderAlias: "coshC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "tanh": {
            isFunction: true,
            shaderAlias: "tanhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "asinh": {
            isFunction: true,
            shaderAlias: "asinhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "acosh": {
            isFunction: true,
            shaderAlias: "acoshC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "atanh": {
            isFunction: true,
            shaderAlias: "atanhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sinp": {
            isFunction: true,
            shaderAlias: "sinpC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "cosp": {
            isFunction: true,
            shaderAlias: "cospC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "re": {
            isFunction: true,
            shaderAlias: "reC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "im": {
            isFunction: true,
            shaderAlias: "imC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "Gamma": {
            isFunction: true,
            shaderAlias: "GammaC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "beta": {
            isFunction: true,
            shaderAlias: "betaC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "min": {
            isFunction: true,
            shaderAlias: "minC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "max": {
            isFunction: true,
            shaderAlias: "maxC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "lerp": {
            isFunction: true,
            shaderAlias: "lerpC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 }, "t": { isFunction: false, type: "complex", index: 2 } },
        },
        "conj": {
            isFunction: true,
            shaderAlias: "conjC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "clamp": {
            isFunction: true,
            shaderAlias: "clampC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 }, "t": { isFunction: false, type: "complex", index: 2 } },
        },
        "frac": {
            isFunction: true,
            shaderAlias: "fracC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "floor": {
            isFunction: true,
            shaderAlias: "floorC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "ceil": {
            isFunction: true,
            shaderAlias: "ceilC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "mod": {
            isFunction: true,
            shaderAlias: "modC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "b": { isFunction: false, type: "complex", index: 0 } },
        },
        "inverseSC": {
            isFunction: true,
            shaderAlias: "inverseSCC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 0 } },
        },
        "sc": {
            isFunction: true,
            shaderAlias: "scC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "planeToP": {
            isFunction: true,
            shaderAlias: "planeToPC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "pToPlane": {
            isFunction: true,
            shaderAlias: "pToPlaneC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "squeeze": {
            isFunction: true,
            shaderAlias: "squeezeC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "coverage": { isFunction: false, type: "complex", index: 1 }, "length": { isFunction: false, type: "complex", index: 2 } },
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
    "powlog": Complex.powlog,
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
    "asinh": Complex.asinh,
    "acosh": Complex.acosh,
    "atanh": Complex.atanh,
    "sinp": Complex.sinp,
    "cosp": Complex.cosp,
    "re": (z) => complex(z.re, 0),
    "im": (z) => complex(z.im, 0),
    "Gamma": Complex.gamma,
    "beta": Complex.beta,
    "min": Complex.min,
    "max": Complex.max,
    "lerp": (z1, z2, t) => Complex.mult(complex(1, 0).sub(t), z1).add(z2.mult(t)),
    "conj": (z) => z.conj(),
    "clamp": (z, min, max) => Complex.clamp(z, min.norm(), max.norm()),
    "frac": Complex.frac,
    "floor": Complex.floor,
    "ceil": Complex.ceil,
    "mod": Complex.mod,
    "inverseSC": (z, p) => complex(1, 0), // not implemented for p5 mode
    "sc": (z, p) => complex(1, 0), // not implemented for p5 mode
    "planeToP": (z, p) => complex(1, 0), // not implemented for p5 mode
    "pToPlane": (z, p) => complex(1, 0), // not implemented for p5 mode
    "squeeze": (z, coverage, length) => Complex.squeeze(z, coverage.re, length.re),
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
    "realBounds": complex(0, 0),
    "imagBounds": complex(0, 0),
};

let valueScope = {};

module.exports = {
    scope, defaultValueScope, valueScope,
};