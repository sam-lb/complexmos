


const scopes = require("../app/scope.js");
const scope = scopes.scope.builtin;

const descriptions = {
    "Gamma": "The Gamma function",
};


const generateDescriptions = () => {
    const keys = Object.keys(scope).sort();
    let result = "";

    for (const key of keys) {
        const description = descriptions[key] ?? "No description given";
        let arguments = "";
        if (scope[key].isFunction) {
            arguments = "Arguments: ";
            const locals = Object.keys(scope[key].locals).sort(local => local.index);
            if (locals.length === 0) arguments += "none";
            for (let i=0; i<locals.length; i++) {
                const local = locals[i];
                arguments += local + ` (${scope[key].locals[local].type})`;
                if (i !== locals.length - 1) arguments += ", ";
            }
            arguments += "<br>";
        }
        result += `<div class="description-entry"><span style="font-weight: bold">${key}</span><br>${arguments}Description:${description}</div>`;
    }

    return result;
};

const outputDescriptions = (targetId) => {
    const targetElement = document.querySelector(`#${targetId}`);
    targetElement.innerHTML = generateDescriptions();
}

const toggleDescriptions = () => {
    const descDiv = document.querySelector("#description-container");
    const collapseBtn = document.querySelector("#collapse-btn");
    if (descDiv.style.display === "none") {
        descDiv.style.display = "block";
        collapseBtn.innerText = "Collapse";
    } else {
        descDiv.style.display = "none";
        collapseBtn.innerText = "Expand";
    }
}

outputDescriptions("description-container");
window.toggleDescriptions = toggleDescriptions;