const scope = {
    builtin: {
        "x": { isFunction: false, },
        "y": { isFunction: false, },
        "xy": { isFunction: false, },
        "x1": { isFunction: false, },
        /* datatypes (for type annotation) */
        "complex": { isFunction: false, },
        "array": { isFunction: false, },
        "matrix": { isFunction: false, },
        "function": { isFunction: false, },
    },
    userGlobal: {
        "f": { isFunction: true, },
    },
}