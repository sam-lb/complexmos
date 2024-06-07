const scope = {
    builtin: {
        "z": { isFunction: false, },
        "sin": { isFunction: true, },
        "i": { isFunction: false, },

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
    "sin": Complex.sin,
    "i": complex(0, 1),
};