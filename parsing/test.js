document.querySelector("#minput-field").focus();

const tracker = new ErrorTracker("error-output", "Tokenization successful!");

const TEST_STRINGS = [
    "\\left[x+\\frac{1}{x},2\\right]\\left(2+3+y\\right)",
    "\\frac{x}{.1y}",
    "clamp\\left(x,lower,upper\\right)=\\max\\left(lower,\\min\\left(upper,x\\right)\\right)",
];

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
            value: complex(2, 0),
        },
        "y": {
            dataType: dataTypes.number,
            value: complex(1, 0),
        },
        "z": {
            dataType: dataTypes.number,
            value: complex(2, 0),
        },
        "xy": {
            dataType: dataTypes.array,
            value: [],
        },
        "i": {
            dataType: dataTypes.number,
            value: complex(0, 1),
        },
        "clamp": {
            dataType: dataTypes.function,
            value: null,
            arguments: [dataTypes.number, dataTypes.number, dataTypes.number],
        },
        "lower": {
            dataType: dataTypes.number,
            value: 1,
        },
        "upper": {
            dataType: dataTypes.number,
            value: 2,
        },
        "max": {
            dataType: dataTypes.function,
            value: null,
            arguments: [dataTypes.number, dataTypes.number],
        },
        "min": {
            dataType: dataTypes.function,
            value: null,
            arguments: [dataTypes.number, dataTypes.number],
        },
    }
};


const LOAD_DICTIONARY = false;
let dictionary;
if (LOAD_DICTIONARY) {
    fetch("https://raw.githubusercontent.com/dwyl/english-words/master/words.txt")
    .then(resp => resp.text())
    .then(text => {
        dictionary = new Trie(text.split("\n"));
        document.querySelector("#prefix-search").disabled = false;
    });
}

function prefixSearch() {
    if (!LOAD_DICTIONARY) return;
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

    if (tokens !== null) {
        let tokenString = "";
        for (let token of tokens.tokens) {
            tokenString += token.toString() + ", ";
        }
        console.log(tokenString.slice(0, -2));

        tracker.clear();
    }
}