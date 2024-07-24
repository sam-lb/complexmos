


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
        result += `<div class="description-entry"><span style="font-weight: bold">${key}</span>: ${description}</div>`;
    }

    return result;
};

const outputDescriptions = (targetId) => {
    const targetElement = document.querySelector(`#${targetId}`);
    targetElement.innerHTML = generateDescriptions();
}

outputDescriptions("description-container");