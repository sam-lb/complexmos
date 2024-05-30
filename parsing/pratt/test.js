

document.querySelector("#minput-field").focus();



const scope = {
    builtin: {
        "x": { isFunction: false, },
        "y": { isFunction: false, },
        "xy": { isFunction: false, },
        "x1": { isFunction: false, },
    },
    userGlobal: {
        "f": { isFunction: true, },
    },
}


function handleKeydown() {
    if (window.event.keyCode === 13) {
        handleSubmit();
    }
}


function handleSubmit() {
    const latex = document.querySelector("#minput-field").value;
    const text = cleanLatex(latex);

    tracker.setTarget("error-output");
    tracker.clear();

    const lexer = new Lexer(text, false);
    lexer.setScope(scope);
    
    lexer.tokenize();
    const tokens = lexer.getTokens();

    if (tracker.hasError) {
        console.error("error occurred during tokenization");
    } else {
        tracker.clear();
        console.log("Tokens: ", tokens);

        const parser = new ExpressionParser(tokens);
        const result = parser.parseExpression();

        if (tracker.hasError) {
            console.error("error occurred during parsing");
        } else {
            console.log("AST:\n\n", result, "\n\n"+result.toString());
        }
    }
}