const { GRADIENTS } = require("./coloring.js");



const buildMenuHTML = (target, id, index, error=null) => {
    const visibility = (error) ? ["none", "flex"] : ["flex", "none"];
    target.innerHTML = `<div id="valid-index-display-${id}" class="expr-index">${index}</div>
                        <div class="side-menu-image" id="icon-container-${id}">
                            <div style="display:${visibility[0]};border-radius:5px;" id="icon-valid-container-${id}">
                                <img src="../data/settings_transparent.png" style="width:25px;height:25px;" onclick="displayOverlayMenu(${id});"></img>
                            </div>
                            <div style="display:${visibility[1]};border-radius:5px;" id="icon-invalid-container-${id}">
                                <img src="../data/error_transparent.png" style="width:25px;height:25px;" onclick="displayOverlayMenu(${id});"></img>
                            </div>
                        </div>`;

};

const modifyMenuHTML = (id, error=null) => {
    const validEl = document.querySelector(`#icon-valid-container-${id}`);
    const invalidEl = document.querySelector(`#icon-invalid-container-${id}`);
    if (error) {
        validEl.style.display = "none";
        document.querySelector(`#icon-container-${id}`).setAttribute("title", error);
        invalidEl.style.display = "flex";
    } else {
        validEl.style.display = "flex";
        invalidEl.style.display = "none";
        document.querySelector(`#icon-container-${id}`).setAttribute("title", "Settings");
    }
};

const modifyExpressionIndex = (id, index) => {
    document.querySelector(`#valid-index-display-${id}`).innerHTML = index.toString();
    fields[id].index = index;
};

function displayOverlayMenu(id) {
    const overlay = document.querySelector("#overlay-menu-container");
    overlay.style.display = "block";
    const additionalSettings = fields[id]["settingsHTML"] ?? "";
    overlay.innerHTML = `
    Settings for expression ${fields[id].index}
    <hr>${additionalSettings}
    `;
}

function hideOverlayMenu() {
    const overlay = document.querySelector("#overlay-menu-container");
    overlay.style.display = "none";
}

function generateSettingsHTML(id) {
    let checked, colorMode, currentGradient;
    if (fields[id]["settingsHTML"] || document.querySelector(`#display-checkbox-${id}`)) {
        checked = fields[id]["displaySettings"]["display"];
        colorMode = fields[id]["displaySettings"]["colorMode"];
        currentGradient = fields[id]["displaySettings"]["gradient"] ?? "monokai";
    } else {
        checked = true;
        colorMode = "default";
        currentGradient = "monokai";
    }
    const checkedString = checked ? " checked" : "";
    fields[id]["displaySettings"] = {
        "display": checked,
        "colorMode": colorMode,
        "gradient": currentGradient,
    };

    let gradientDropdown = "";
    if (["gradient", "gradient-discrete"].includes(colorMode)) {
        let gradientOptions = "";
        for (let gradient of Object.keys(GRADIENTS).sort()) {
            const checkedGrad = (gradient === currentGradient) ? " selected" : "";
            gradientOptions += `<option value="${gradient}" ${checkedGrad}>${gradient}</option>`;
        }
        gradientDropdown = `<label for="gradient-dropdown-${id}">Gradient</label>
                            <select id="gradient-dropdown-${id}" onchange="plot.setGradientMode(${id}, true);">
                                ${gradientOptions}
                            </select><br>`;
    }

    return `<label for="display-checkbox-${id}">Display?</label>
    <input type="checkbox" id="display-checkbox-${id}" oninput="plot.toggleDisplay(${id});"${checkedString}><br>
    <label for="display-coloring-dropdown-${id}">Coloring mode</label>
    <select id="display-coloring-dropdown-${id}" onchange="plot.setColorMode(${id});">
        <option value="default">Default (HSV rainbow)</option>
        <option value="default-discrete">Discrete Default</option>
        <option value="gradient">Gradient</option>
        <option value="gradient-discrete">Discrete Gradient</option>
        <option value="image-repeat">Image (repeated)</option>
        <option value="image-stretch">Image (stretch)</option>
    </select><br>
    ${gradientDropdown}
    <div style="display:flex;flex-direction:row;">
        <button class="sick-btn-design" onclick="plot.uploadImage(${id});">Upload image</button>
        <div>${fields[id]["imageFile"] ?? "none selected"}</div>
    </div>
    `.replace(`value="${colorMode}">`, `value="${colorMode}" selected>`); // ah yes
}

function handleSlider(id) {
    const slider = document.querySelector(`#slider-${id}`);
    const variable = fields[id].field;
    const assignment = variable.latex().split("=")[0];
    variable.latex(`${assignment}=${parseFloat(slider.value)}`);
}

function bottomHTML(target, bounds, id) {
    const div = document.querySelector(`#${target}`);
    const oldContainer = document.querySelector(`#slider-container-${id}`);
    if (bounds === null) {
        if (sliderFields[id]) delete sliderFields[id];
        if (oldContainer) div.removeChild(oldContainer);
        return;
    }
    if (oldContainer) {
        const slider = document.querySelector(`#slider-${id}`);
        const calculatedBounds = (sliderFields[id]["getBounds"] ?? (() => [0, 1]))();
        slider.setAttribute("min", calculatedBounds[0].toString());
        slider.setAttribute("max", calculatedBounds[1].toString());
        slider.setAttribute("step", `${(calculatedBounds[1] - calculatedBounds[0]) / 100}`);
        slider.value = fields[id].field.latex().split("=")[1];
        return;
    }

    const container = document.createElement("div");
    container.setAttribute("class", "slider-container");
    container.setAttribute("id", `slider-container-${id}`);
    const startSpan = document.createElement("span");
    
    const slider = document.createElement("input");
    slider.setAttribute("type", "range");
    slider.setAttribute("min", `${bounds.min}`);
    slider.setAttribute("max", `${bounds.max}`);
    const step = (bounds.max - bounds.min) / 100;
    slider.setAttribute("step", `${step}`);
    slider.setAttribute("id", `slider-${id}`);
    slider.setAttribute("class", "variable-slider");
    slider.setAttribute("value", fields[id].field.latex().split("=")[1]);
    slider.setAttribute("oninput", `handleSlider(${id});`);

    const endSpan = document.createElement("span");
    startSpan.setAttribute("id", `start-field-${id}`);
    endSpan.setAttribute("id", `end-field-${id}`);

    container.appendChild(startSpan);
    container.appendChild(slider);
    container.appendChild(endSpan);
    div.appendChild(container);

    const startField = MQ.MathField(startSpan, {});
    const endField = MQ.MathField(endSpan, {});

    sliderFields[id] = {
        "min": startField,
        "max": endField,
    };

    startField.latex(`${bounds.min}`);
    endField.latex(`${bounds.max}`);
}

function addField(parent=null) {
    /** add new math input field. parent: parent element */

    const newIndex = Object.keys(fields).length + 1;

    const newDiv = document.createElement("div");
    newDiv.setAttribute("class", "math-input-div-container");

    const subDiv = document.createElement("div");
    subDiv.setAttribute("class", "math-input-div");

    const newSpan = document.createElement("span");
    newSpan.setAttribute("class", "math-input");

    const newField = MQ.MathField(newSpan, {});
    newDiv.setAttribute("id", `math-input-div-container-${newField.id}`);

    const newMenu = document.createElement("div");
    newMenu.setAttribute("class", "math-input-side-menu");
    newMenu.setAttribute("id", `math-input-side-menu-${newField.id}`);
    buildMenuHTML(newMenu, newField.id, newIndex, null);

    const bottomDiv = document.createElement("div");
    bottomDiv.setAttribute("id", `math-input-bottom-div-${newField.id}`);
    subDiv.appendChild(newMenu);
    subDiv.appendChild(newSpan);
    newDiv.appendChild(subDiv);
    newDiv.appendChild(bottomDiv);

    if (parent === null) {
        const container = document.querySelector("#math-input-container");
        container.appendChild(newDiv);

        fields[newField.id] = {
            id: newField.id,
            field: newField,
            last: null,
            next: null,
            container: newDiv,
            displaySettings: {},
            index: newIndex,
        };
    } else {
        const lastDiv = document.querySelector(`#math-input-div-container-${parent.id}`);
        lastDiv.after(newDiv);

        fields[newField.id] = {
            id: newField.id,
            field: newField,
            last: parent.field,
            next: parent.next,
            container: newDiv,
            displaySettings: {},
            index: newIndex,
        };
        fields[parent.field.id].next = newField;

        advance(parent.field.id, 1);
    }

    return newField.id;
}

function deleteField(id, preserve=true) {
    if (preserve && Object.keys(fields).length === 1) return; // at least one field has to remain
    const entry = fields[id];
    if (entry.next !== null) {
        if (entry.last !== null) {
            fields[entry.next.id]["last"] = entry.last;
            fields[entry.last.id]["next"] = entry.next;
        } else {
        fields[entry.next.id]["last"] = null;
        }
    } else {
        if (entry.last !== null) {
            fields[entry.last.id]["next"] = null;
        } else {
            // there are no fields left
        }
    }
    if (preserve) advance(id, (entry.last === null) ? 1 : -1);

    entry.container.parentNode.removeChild(entry.container);
    const index = fields[id].index;
    delete fields[id];
    for (const exprID in fields) {
        if (fields[exprID].index > index) {
            modifyExpressionIndex(exprID, fields[exprID].index - 1);
        }        
    }
}

function advance(id, direction) {
    const entry = fields[id];
    if (direction === -1 && entry.last !== null) {
        entry.last.focus();
        entry.last.moveToRightEnd();
    } else if (direction === 1) {
        if (entry.next !== null) {
            entry.next.focus();
            entry.next.moveToRightEnd();
        } else {
            addField(entry);
        }
    }
}

function highlightExpression(id) {
    for (const exprID in fields) {
        const icon = document.querySelector(`#icon-valid-container-${exprID}`) || document.querySelector(`#icon-invalid-container-${exprID}`);
        if (id === exprID) {
            icon.style.outline = "thick double darkslategray";
        } else {
            icon.style.outline = "none";
        }
    }
}

function tabSwitch(tab) {
    const plane = document.querySelector("#ui-header-plane");
    const sphere = document.querySelector("#ui-header-sphere");
    const cube = document.querySelector("#ui-header-cube");
    if (tab === 0) {
        plane.style.backgroundColor = "white";
        sphere.style.backgroundColor = "lightgray";
        cube.style.backgroundColor = "lightgray";
        plot.setMode(Plot.modes.PLANE);
    } else if (tab === 1) {
        plane.style.backgroundColor = "lightgray";
        sphere.style.backgroundColor = "white";
        cube.style.backgroundColor = "lightgray";
        plot.setMode(Plot.modes.SPHERE);
    } else {
        plane.style.backgroundColor = "lightgray";
        sphere.style.backgroundColor = "lightgray";
        cube.style.backgroundColor = "white";
        plot.setMode(Plot.modes.CUBE);
    }
}

function toggleSettingsPopup() {
    const popup = document.querySelector("#settings-popup");
    if (popup.style.display !== "flex") {
        popup.style.display = "flex";
    } else {
        popup.style.display = "none";
    }
}



module.exports = {
    modifyMenuHTML, displayOverlayMenu, hideOverlayMenu, generateSettingsHTML,
    handleSlider, bottomHTML, addField, deleteField,
    advance, highlightExpression, tabSwitch, toggleSettingsPopup,
};