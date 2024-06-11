const scope = {
    builtin: {
        "z": { isFunction: false, },

        "norm": { isFunction: true, },
        "normSq": { isFunction: true, },
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

        "i": { isFunction: false, },
        "pi": { isFunction: false, },
        "tau": { isFunction: false, },
        "e": { isFunction: false, },

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
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
};