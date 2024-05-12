document.querySelector("#minput-field").focus();

const tracker = new ErrorTracker("error-output", "Tokenization successful!");

const scope = {
    "builtin": {
        "sin": {
            dataType: dataTypes.function,
            value: Complex.sin,
            arguments: [dataTypes.number],
        },
        "sinh": {
            dataType: dataTypes.function,
            value: Complex.sinh,
            arguments: [dataTypes.number],
        },
        "matmul": {
            dataType: dataTypes.function,
            value: Matrix.multiply,
            arguments: [dataTypes.matrix, dataTypes.matrix],
        },
        "x": {
            dataType: dataTypes.number,
            value: 0,
        },
        "y": {
            dataType: dataTypes.number,
            value: 1,
        },
        "z": {
            dataType: dataTypes.number,
            value: 2,
        },
        "xy": {
            dataType: dataTypes.array,
            value: [],
        },
    }
};


function handleKeydown() {
    if (window.event.keyCode === 13) {
        handleSubmit();
    }
}


function handleSubmit() {
    const latex = document.querySelector("#minput-field").value;
    const text = cleanLatex(latex);
    const tokens = tokenize(text, tracker, scope);
    console.log(tokens?.tokens);

    if (tokens !== null) {
        tracker.clear();
    }
}