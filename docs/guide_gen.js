


const scopes = require("../app/scope.js");
const scope = scopes.scope.builtin;

const descriptions = {

};


const generateDescriptions = () => {
    const keys = Object.keys(scope).sort();
    let result = "";

    for (const key of keys) {
        const description = descriptions[key] ?? "No description given";
        result += `<div><span style="font-weight: bold">${key}</span>: ${description}</div>`;
    }

    return result;
};

module.exports = {
    generateDescriptions
};