const { tracker } = require("../parsing/errors.js");
const { scope, defaultValueScope, valueScope } = require("./scope.js");
const { evaluate } = require("./evaluator.js");
const { classifyInput, classifySliderInput, validateLines, populateUserScope, validateAST } = require("./expression_processor.js");
const { translateToGLSL, colorGLSLFromSettings } = require("./translator.js");
const { VariableDefinition, FunctionDefinition } = require("./input_expressions.js");
const {
    menuHTML, displayOverlayMenu, generateSettingsHTML,
    handleSlider, bottomHTML, addField, deleteField, advance,
    highlightExpression, tabSwitch, toggleSettingsPopup,
} = require("./domhandler.js");
const { DomainColoring } = require("./plottable.js");
const { loadShaders, preload } = require("./data_loading.js");
const { registerMouseEvents, windowResized } = require("./event_handling.js");
const { Plot } = require("./plot.js");

p5.disableFriendlyErrors = true;
window.plot = undefined;
window.RENDERER = "WebGL";


/** MathQuill handling */


window.MQ = MathQuill.getInterface(2);
const opsString = Object.keys(scope.builtin).filter(x => x.length > 1).join(" ");
window.fields = {};
window.sliderFields = {};

function handleDisplayToggles(lines) {
    if (!lines) return;
    const plottableIDs = [];
    for (const line of lines) {
        if (!(line instanceof FunctionDefinition)) {
            fields[line.id]["settingsHTML"] = "";
            fields[line.id]["displaySettings"]["display"] = false;
            continue;
        }
        const locals = Object.keys(scope.userGlobal[line.name].locals);
        if (locals.length !== 1 || locals[0] !== "z") {
            fields[line.id]["settingsHTML"] = "";
            fields[line.id]["displaySettings"]["display"] = false;
            continue;
        }

        plottableIDs.push(line.id);
    }



    if (plottableIDs.length === 0) return;
    for (id of plottableIDs) {
        fields[id]["settingsHTML"] = generateSettingsHTML(id);;
    }
}


function debounceWrapper(func, interval, initialTimer) {
    let timer = initialTimer;
    return function() {
        const context = this;
        const args = arguments;
        clearTimeout(timer);
        timer = setTimeout(() => {func.apply(context, args);}, interval);
    };
}

window.getCallbacks = (id) => {
    const sideMenu = document.querySelector(`#math-input-side-menu-${id}`);

    const callback = (message, target) => {
        sideMenu.innerHTML = menuHTML(id, message);
    };

    const successCallback = (message, target) => {
        sideMenu.innerHTML = menuHTML(id);        
    };

    return {
        callback: callback,
        successCallback: successCallback,
    };
}

function validateInput() {
    populateUserScope(fields);
    if (tracker.hasError) return null;

    const lines = classifyInput(fields);
    if (tracker.hasError) return null;

    validateLines(lines);
    if (tracker.hasError) return null;

    for (const line of Array.prototype.concat(lines["functions"], lines["variables"], lines["evaluatables"])) {
        const callbacks = getCallbacks(line.id);
        tracker.setCallback(callbacks.callback);
        tracker.setSuccessCallback(callbacks.successCallback);
        line.buildAST();
        if (tracker.hasError) return null;
        validateAST(line.ast);
        if (tracker.hasError) return null;
    }

    const varsAndFuncs = lines["functions"].concat(lines["variables"]);
    const noNumbers = (string) => !(["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"].some((x) => string.includes(x)));
    let newOpsString = opsString + " " + varsAndFuncs.filter(line => line.name.length > 1 && noNumbers(line.name)).map(line => line.name).join(" ");
    if (newOpsString[newOpsString.length-1] === " ") newOpsString = newOpsString.slice(0, -1);
    MQ.config({
        autoOperatorNames: newOpsString,
    });

    return varsAndFuncs;
}

function validateSliderInput() {
    const sliderLines = classifySliderInput(sliderFields);
    if (tracker.hasError) return null;

    for (const line of sliderLines) {
        const minLine = line[0];
        const maxLine = line[1];

        const callbacks = getCallbacks(minLine.id);
        tracker.setCallback(callbacks.callback);
        tracker.setSuccessCallback(callbacks.successCallback);

        minLine.buildAST();
        if (tracker.hasError) return null;
        validateAST(minLine.ast);
        if (tracker.hasError) return null;

        maxLine.buildAST();
        if (tracker.hasError) return null;
        validateAST(maxLine.ast);
        if (tracker.hasError) return null;
    }

    return sliderLines;
}

function setBoundCalculators(sliderLines) {
    if (sliderLines === null) return;
    for (const line of sliderLines) {
        const minLine = line[0], maxLine = line[1];
        const minEval = evaluate(minLine.ast), maxEval = evaluate(maxLine.ast);
        sliderFields[minLine.id]["getBounds"] = () => {
            const minBound = minEval.call();
            const maxBound = maxEval.call();
            return (minBound && maxBound) ? [minEval.call().re, maxEval.call().re] : [0, 1];
        };
    }
}

function populateValueScope(lines) {
    Object.keys(valueScope).forEach(key => delete valueScope[key]);
    Object.keys(defaultValueScope).forEach(key => valueScope[key] = defaultValueScope[key]);
    if (lines === null) return;
    for (const line of lines) {
        if (scope.userGlobal[line.name].isFunction) {
            const locals = scope.userGlobal[line.name].locals;
            valueScope[line.name] = evaluate(line.ast, Object.keys(locals).sort(key => locals[key].index));
        } else {
            valueScope[line.name] = evaluate(line.ast);
        }
    }
}

function configureRenderers(lines) {
    if (RENDERER === "WebGL") {
        if (lines === null) {
            plot.setShaderReplacement(null);
            plot.setDisplayReplacement(null);
        } else {
            const emittedGLSL = translateToGLSL(lines.slice());
            if (emittedGLSL) {
                plot.setShaderReplacement(emittedGLSL);
                const displayName = pickDisplay(lines);
                highlightExpression(displayName?.id);
                if (displayName) {
                    plot.setDisplayReplacement(displayName.name, colorGLSLFromSettings(displayName.id));
                } else {
                    plot.setDisplayReplacement(null);
                }
            } else {
                plot.setShaderReplacement(null);
                plot.setDisplayReplacement(null);
            }
        }
    } else {
        // use evaluate() and scope.userGlobal to populate valueScope
        plot.clear();
        if (!lines) return;
        let displayName = pickDisplay(lines);
        highlightExpression(displayName?.id);
        if (!displayName) return;
        displayName = displayName.name;
        plot.addPlottable(new DomainColoring((z) => valueScope[displayName].call({z:z}),));
    }
}

function addSliders(lines) {
    if (lines === null) return;
    for (const line of lines) {
        if (line instanceof VariableDefinition) {
            const bounds = line.sliderBounds(fields);
            bottomHTML(`math-input-bottom-div-${line.id}`, bounds, line.id);
        } else {
            bottomHTML(`math-input-bottom-div-${line.id}`, null, line.id);
        }
    }
}

function pickDisplay(lines) {
    const rev = Object.keys(fields).sort((a, b) => parseInt(b) - parseInt(a));
    for (const id of rev) {
        if (fields[id]["displaySettings"]["display"]) {
            return {id: id, name: lines.filter(l => l.id === id)[0]?.name };
        }
    }
}

function fieldEditHandler(mathField) {
    if (mathField === null || fields[mathField.id]) {
        const lines = validateInput();
        const sliderLines = validateSliderInput();
        handleDisplayToggles(lines);
        populateValueScope(lines);
        setBoundCalculators(sliderLines);
        addSliders(lines);
        configureRenderers(lines);
    } else {
        // slider field
        fieldEditHandler(null);
    }
}

const firstField = addField();
fields[firstField].field.focus();

MQ.config({
    autoCommands: "pi sqrt tau alpha beta Gamma",
    supSubsRequireOperand: true,
    charsThatBreakOutOfSupSub: "",
    autoOperatorNames: opsString,
    handlers: {
        moveOutOf: (direction, mathField) => {
            if (!(direction === MQ.L || direction === MQ.R)) mathField.moveToLeftEnd();
        },
        enter: (mathField) => {
            mathField.moveToLeftEnd();
            advance(mathField.id, 1);
        },
        downOutOf: (mathField) => { 
            mathField.moveToLeftEnd();
            advance(mathField.id, 1);
        },
        upOutOf: (mathField) => {
            mathField.moveToLeftEnd();
            advance(mathField.id, -1);
        },
        deleteOutOf: (direction, mathField) => { if (direction === MQ.L) deleteField(mathField.id); },
        edit: (mathField) => {
            if (fields[mathField.id]) {
                debounceWrapper(fieldEditHandler, 500, -1)(mathField);
            } else {
                fieldEditHandler(mathField);
            }
        }
    },
});


/** ******************************** */

function setupP5(offset=null) {
    const canvasDiv = document.querySelector("#canvas-div");
    canvasDiv.innerHTML = "";
    const canvas = createCanvas(canvasDiv.offsetWidth, canvasDiv.offsetHeight);
    canvas.parent("canvas-div");

    const mode = plot?.mode ?? Plot.modes.PLANE;
    plot = new Plot(canvasDiv.offsetWidth, canvasDiv.offsetHeight, null, mode, false);
    if (offset) plot.pan(offset)
    tabSwitch(plot.mode-1);
    fieldEditHandler(null);
}

function reglLoaded(err, regl, shaders, offset) {
    if (err) {
        console.warn(`Could not load WebGL! Maybe your browser doesn't support it? Using vanilla canvas instead. Specific error: ${err}`);
        RENDERER = "p5";
        setupP5();
        fieldEditHandler(null);
        return;
    }

    console.log("regl loaded!");

    regl.on("lost", contextLost);
    const canvasDiv = document.querySelector("#canvas-div");
    const mode = plot?.mode ?? Plot.modes.PLANE;
    plot = new Plot(canvasDiv.offsetWidth, canvasDiv.offsetHeight, null, mode, false, regl, shaders);
    if (offset) plot.pan(offset)
    tabSwitch(plot.mode-1);    
    fieldEditHandler(null);
}

function shadersLoaded(shaders, offset=null) {
    if (plot?.reglInstance) {
        plot.reglInstance._refresh();
        reglLoaded(null, plot.reglInstance, shaders, offset);
    } else {
        // remove the loading shaders message
        document.querySelector("#canvas-div").innerHTML = "";
        require("regl")({
            container: "#canvas-div",
            onDone: (err, regl) => reglLoaded(err, regl, shaders, offset),
        });
    }
}

function setupWebGL(offset=null) {
    if (plot?.shaders) {
        shadersLoaded(plot.shaders, offset);
    } else {
        loadShaders().then(shaders => shadersLoaded(shaders, offset));
    }
}

function setup(offset=null) {
    if (RENDERER === "p5") {
        setupP5(offset);
    } else {
        setupWebGL(offset); 
    }

    registerMouseEvents();
}

function draw() {
    if (plot) plot.update();
}

function setRenderer() {
    const renderer = document.querySelector("#webgl-toggle").checked;
    if (RENDERER === "WebGL" && !renderer) {
        RENDERER = "p5";
        setupP5();
    } else if (RENDERER === "p5" && renderer) {
        RENDERER = "WebGL";
        setupWebGL();
    }
}

function contextLost() {
    const webGLToggle = document.querySelector("#webgl-toggle");
    webGLToggle.checked = false;
    webGLToggle.disabled = true;
    alert("Your WebGL context has been lost :(");

    setRenderer();
}

// there might be a better way to do this, but it's actually fine
window.resizeDebounced = debounceWrapper(() => windowResized(setup), 250, -1);
window.preload = preload;
window.displayOverlayMenu = displayOverlayMenu;
window.tabSwitch = tabSwitch;
window.draw = draw;
window.addEventListener("resize", resizeDebounced);
window.toggleSettingsPopup = toggleSettingsPopup;
window.setRenderer = setRenderer;
window.handleSlider = handleSlider;
window.fieldEditHandler = fieldEditHandler;
window.Plot = Plot;

window.onload = () => {
    const aspect = window.innerWidth / window.innerHeight;
    if (aspect <= 1) {
        alert("Rotate your device for the best experience. Note that features may not work as expected on mobile.");
    } else if (aspect < 4 / 3) {
        alert(`It is recommended to use this app on a device with 4:3 (1.33:1) aspect ratio or greater. Your device aspect ratio: ${Math.floor(100 * aspect) / 100}:1`);
    }

    fetch("../package.json").then(response => response.json().then(json => {
        document.querySelector("#version-div").innerHTML = `<a href="https://github.com/sam-lb/complexmos/releases" target="_blank">v${json.version}</a>`;
    }));

    setup();
};