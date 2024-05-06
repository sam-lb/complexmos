document.querySelector("#minput-field").focus();

const tracker = new ErrorTracker("error-output", "Tokenization successful!");

const scope = {
    "builtin": [
        "sin", "sinh", "max", "x"
    ],
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
    console.log(tokens);

    if (tokens !== null) {
        tracker.clear();
    }
}