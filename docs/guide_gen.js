


const scopes = require("../app/scope.js");
const scope = scopes.scope.builtin;

const descriptions = {
    "Gamma": "The Gamma function \\( \\Gamma(z) \\).",

    "acos": "The inverse cosine function.",

    "arg": "The argument or angle of a complex number.",

    "array": "The array datatype. This is currently not yet implemented, but the keyword is already reserved in the calculator.",

    "asin": "The inverse sine function.",

    "atan": "The inverse tangent function.",

    "atanh": "The inverse hyperbolic tangent function.",

    "beta": "The beta function \\( \\beta(z, w) \\).",

    "clamp": "Make sure that \\( |w| \\leq |z| \\leq |t|  \\). If so, return \\( z \\). If not, return \\( z \\) scaled to the norm of the closer bound. For example, \\( \\operatorname{clamp}(1.1, 0, 1)=1 \\).",

    "complex": "The complex datatype. It represents a single number of the form \\( a+ib \\).",

    "conj": "The complex conjugate of a complex number, that is, the reflection of the number across the real axis.",

    "cos": "The cosine function.",

    "cosh": "The hyperbolic cosine function.",

    "e": "Euler's constant, the base of the natural logarithm, \\( e \\approx 2.718 \\).",

    "exp": "The exponential function.",

    "frac": "Returns the fractional part (the part after the decimal point) of both the real and imaginary components of a complex number. For example, \\( \\operatorname{frac}(1.123 + 4.567i)=0.123+0.567i \\).",

    "function": "The function datatype",
    
    "i": "The imaginary unit \\( i := x \\in \\mathbb{R}\\left[x\\right]/(x^2+1) \\).",

    "im": "The imaginary part of a complex number.",

    "imagBounds": "The current y-axis bounds in plane mode. The format is \\( \\operatorname{imagBounds} = \\operatorname{yMin} + \\operatorname{yMax}\\cdot i \\).",

    "inv": "The multiplicative inverse of a complex number.",

    "inverseSC": "Construction of the inverse Schwarz-Christoffel map for a regular \\( p \\)-gon. This conformally maps the unit disk to the \\( p \\)-gon.",

    "lerp": "Linearly interpolates \\( z \\) between \\( w \\) and \\( t \\).",

    "ln": "The natural (base \\( e \\) logarithm.",

    "matrix": "The matrix datatype. This is currently not yet implemented, but the keyword is already reserved in the calculator.",
    
    "max": "Returns the argument with the greater norm.",

    "min": "Returns the argument with the lesser norm.",
    
    "norm": "The Euclidean norm (magnitude) of a complex number.",

    "normSq": "The square of the Euclidean norm (magnitude) of a complex number.",

    "pToPlane": "A conformal map from a regular unit-radius \\( p \\)-gon to the entirety of the complex plane.",

    "pi": "The mathematical half-circle constant \\( \\pi \\approx 3.142 \\).",

    "planeToP": "A conformal map from the complex plane to a regular unit-radius \\( p \\)-gon.",

    "re": "The real part of a complex number.",

    "realBounds": "The current x-axis bounds in plane mode. The format is \\( \\operatorname{imagBounds} = \\operatorname{xMin} + \\operatorname{xMax}\\cdot i \\).",

    "sc": "Construction of the Schwarz-Christoffel map for a regular \\( p \\)-gon. This conformlly maps a regular \\( p \\)-gon to the unit disk.",

    "sin": "The sine function.",

    "sinh": "The hyperbolic sine function.",

    "sqrt": "The square root function.",

    "squeeze": `Smoothly compresses a \\( 2\\cdot\\operatorname{length} \\) by \\( 2\\cdot\\operatorname{length} \\) square centered at the origin to the unit circle. The incircle of the square is mapped to a circle of radius \\( 0 \\leq \\operatorname{coverage} \\leq 1 \\). You can get a feel for what this does by trying it out on the identity function. See <a href="https://www.desmos.com/calculator/imksvk6d42">here</a> for implementation details.`,

    "tan": "The tangent function.",

    "tanh": "The hyperbolic tangent function.",

    "tau": "The mathematical circle constant \\( \\tau \\approx 6.283 \\).",

    "z": "A parameter representing a complex number that ranges across the entire plane.",
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
        result += `<div class="description-entry"><span style="font-weight: bold">${key}</span><br>${arguments}Description: ${description}</div>`;
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