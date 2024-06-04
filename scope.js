const scope = {
    builtin: {
        "z": { isFunction: false, },

        /* datatypes (for type annotation) */
        "complex": { isFunction: false, isType: true, },
        "array": { isFunction: false, isType: true, },
        "matrix": { isFunction: false, isType: true, },
        "function": { isFunction: false, isType: true, },
    },

    userGlobal: {
    },
}