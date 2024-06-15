const scope = {
    builtin: {
        "z": { isFunction: false, },

        "norm": { isFunction: true, },
        "normSq": { isFunction: true, },
        "arg": { isFunction: true, },
        "inv": { isFunction: true, },
        "exp": { isFunction: true, },
        "ln": { isFunction: true, },
        "sqrt": { isFunction: true, },
        "sin": { isFunction: true, },
        "cos": { isFunction: true, },
        "tan": { isFunction: true, },
        "sinh": { isFunction: true, },
        "cosh": { isFunction: true, },
        "tanh": { isFunction: true, },
        "re": { isFunction: true, },
        "im": { isFunction: true, },
        "Gamma": { isFunction: true, },
        "beta": { isFunction: true, },
        "min": { isFunction: true, },
        "max": { isFunction: true, },

        "i": { isFunction: false, },
        "pi": { isFunction: false, },
        "tau": { isFunction: false, },
        "e": { isFunction: false, },
        "realBounds": { isFunction: false, },
        "imagBounds": { isFunction: false, },

        /* datatypes (for type annotation) */
        "complex": { isFunction: false, isType: true, },
        "array": { isFunction: false, isType: true, },
        "matrix": { isFunction: false, isType: true, },
        "function": { isFunction: false, isType: true, },
    },

    userGlobal: {
    },
}


const valueScope = {
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
    "sinh": Complex.sinh,
    "cosh": Complex.cosh,
    "tanh": Complex.tanh,
    "re": (z) => complex(z.re, 0),
    "im": (z) => complex(z.im, 0),
    "Gamma": Complex.gamma,
    "beta": Complex.beta,
    "min": Complex.min,
    "max": Complex.max,
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
    "realBounds": complex(0, 0),
    "imagBounds": complex(0, 0),
};