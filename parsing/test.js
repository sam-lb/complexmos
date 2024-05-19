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


let dictionary;
fetch("https://raw.githubusercontent.com/dwyl/english-words/master/words.txt")
    .then(resp => resp.text())
    .then(text => {
        dictionary = new Trie(text.split("\n"));
        document.querySelector("#prefix-search").disabled = false;
    });


function prefixSearch() {
    if (window.event.keyCode === 13) {
        const prefix = document.querySelector("#prefix-search").value;
        const results = dictionary.prefixSearch(prefix, 100);

        let htmlString = "";
        const displayCount = Math.min(100, results.length);
        for (let i=0; i<displayCount; i++) {
            const result = results[i];
            htmlString += `<li>${result}</li>`;
        }
        document.querySelector("#results-count").innerHTML = `${results.length} results found for prefix '${prefix}'`;
        document.querySelector("#results-list").innerHTML = htmlString;
    }
}


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