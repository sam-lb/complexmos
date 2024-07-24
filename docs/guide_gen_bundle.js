(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){

const { complex, Complex } = require("../math/complex.js");


const scope = {
    builtin: {
        "z": {
            isFunction: false,
            shaderAlias: "z",
            isParameter: true,
        },
        // "t": {
        //     isFunction: false,
        //     shaderAlias: "t",
        //     isParameter: true,
        // },

        "norm": {
            isFunction: true,
            shaderAlias: "normC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "normSq": {
            isFunction: true,
            shaderAlias: "normSqC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "arg": {
            isFunction: true,
            shaderAlias: "argC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "inv": {
            isFunction: true,
            shaderAlias: "invC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "exp": { 
            isFunction: true,
            shaderAlias: "expC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "ln": { 
            isFunction: true,
            shaderAlias: "lnC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sqrt": { 
            isFunction: true,
            shaderAlias: "sqrtC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sin": { 
            isFunction: true,
            shaderAlias: "sinC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "cos": { 
            isFunction: true,
            shaderAlias: "cosC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "tan": {
            isFunction: true,
            shaderAlias: "tanC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "asin": {
            isFunction: true,
            shaderAlias: "asinC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "acos": {
            isFunction: true,
            shaderAlias: "acosC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "atan": {
            isFunction: true,
            shaderAlias: "atanC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "sinh": {
            isFunction: true,
            shaderAlias: "sinhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "cosh": {
            isFunction: true,
            shaderAlias: "coshC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "tanh": {
            isFunction: true,
            shaderAlias: "tanhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "atanh": {
            isFunction: true,
            shaderAlias: "atanhC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "re": {
            isFunction: true,
            shaderAlias: "reC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "im": {
            isFunction: true,
            shaderAlias: "imC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "Gamma": {
            isFunction: true,
            shaderAlias: "GammaC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "beta": {
            isFunction: true,
            shaderAlias: "betaC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "min": {
            isFunction: true,
            shaderAlias: "minC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "max": {
            isFunction: true,
            shaderAlias: "maxC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 } },
        },
        "lerp": {
            isFunction: true,
            shaderAlias: "lerpC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 }, "t": { isFunction: false, type: "complex", index: 2 } },
        },
        "conj": {
            isFunction: true,
            shaderAlias: "conjC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "clamp": {
            isFunction: true,
            shaderAlias: "clampC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "w": { isFunction: false, type: "complex", index: 1 }, "t": { isFunction: false, type: "complex", index: 2 } },
        },
        "frac": {
            isFunction: true,
            shaderAlias: "fracC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 } },
        },
        "inverseSC": {
            isFunction: true,
            shaderAlias: "inverseSCC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 0 } },
        },
        "sc": {
            isFunction: true,
            shaderAlias: "scC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "planeToP": {
            isFunction: true,
            shaderAlias: "planeToPC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "pToPlane": {
            isFunction: true,
            shaderAlias: "pToPlaneC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "p": { isFunction: false, type: "complex", index: 1 } },
        },
        "squeeze": {
            isFunction: true,
            shaderAlias: "squeezeC",
            locals: { "z": { isFunction: false, type: "complex", index: 0 }, "coverage": { isFunction: false, type: "complex", index: 1 }, "length": { isFunction: false, type: "complex", index: 2 } },
        },

        "i": {
            isFunction: false,
            shaderAlias: "i",
        },
        "pi": {
            isFunction: false,
            shaderAlias: "piC",
        },
        "tau": {
            isFunction: false,
            shaderAlias: "tpiC",
        },
        "e": {
            isFunction: false,
            shaderAlias: "e",
        },
        "realBounds": {
            isFunction: false,
            shaderAlias: "xBounds",
        },
        "imagBounds": {
            isFunction: false,
            shaderAlias: "yBounds",
        },

        /* datatypes (for type annotation) */
        "complex": { isFunction: false, isType: true, },
        "array": { isFunction: false, isType: true, },
        "matrix": { isFunction: false, isType: true, },
        "function": { isFunction: false, isType: true, },
    },

    userGlobal: {
    },
}


const defaultValueScope = {
    "z": complex(1, 0),
    "norm": (z) => complex(z.norm(), 0),
    "normSq": (z) => complex(z.normSq(), 0),
    "arg": (z) => complex(z.arg(), 0),
    "inv": Complex.inv,
    "exp": Complex.exp,
    "ln": Complex.ln,
    "sqrt": Complex.sqrt,
    "sin": Complex.sin,
    "cos": Complex.cos,
    "tan": Complex.tan,
    "acos": Complex.acos,
    "asin": Complex.asin,
    "atan": Complex.atan,
    "sinh": Complex.sinh,
    "cosh": Complex.cosh,
    "tanh": Complex.tanh,
    "atanh": (z) => complex(1, 0), // not implemented for p5 mode
    "re": (z) => complex(z.re, 0),
    "im": (z) => complex(z.im, 0),
    "Gamma": Complex.gamma,
    "beta": Complex.beta,
    "min": Complex.min,
    "max": Complex.max,
    "lerp": (z1, z2, t) => Complex.mult(complex(1, 0).sub(t), z1).add(z2.mult(t)),
    "conj": (z) => z.conj(),
    "clamp": (z, min, max) => Complex.clamp(z, min.norm(), max.norm()),
    "frac": Complex.frac,
    "inverseSC": (z, p) => complex(1, 0), // not implemented for p5 mode
    "sc": (z, p) => complex(1, 0), // not implemented for p5 mode
    "planeToP": (z, p) => complex(1, 0), // not implemented for p5 mode
    "pToPlane": (z, p) => complex(1, 0), // not implemented for p5 mode
    "squeeze": (z, coverage, length) => complex(1, 0), // not implemented for p5 mode
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
    "realBounds": complex(0, 0),
    "imagBounds": complex(0, 0),
};

let valueScope = {};

module.exports = {
    scope, defaultValueScope, valueScope,
};
},{"../math/complex.js":3}],2:[function(require,module,exports){



const scopes = require("../app/scope.js");
const scope = scopes.scope.builtin;

const descriptions = {
    "Gamma": "The Gamma function \\( \\Gamma(z) \\)",
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
},{"../app/scope.js":1}],3:[function(require,module,exports){
const EPSILON = 0.000001;

const pValues = [
	0.99999999999980993,
	676.5203681218851,
	-1259.1392167224028,
	771.32342877765313,
	-176.61502916214059,
	12.507343278686905,
	-0.13857109526572012,
	9.9843695780195716e-6,
	1.5056327351493116e-7
];


const frac = (x) => {
	return x - (x > 0 ? Math.floor(x) : Math.ceil(x));
};

const clamp = (x, min, max) => {
	return Math.min(max, Math.max(min, x));
};


class Complex {

	/*
	Class for representing complex numbers of the form a + bi
	*/

	constructor(real, imaginary) {
		this.re = real;
		this.im = imaginary;
	}

	conj() {
		/* Computes the complex conjugate */
		return new Complex(this.re, -this.im);
	}

	norm() {
		/* Computes the norm (modulus), as a real number */
		return Math.sqrt(this.re * this.re + this.im * this.im);
	}

	normSq() {
		/* Computes the square of the norm (modulus), as a real number */
		return this.re * this.re + this.im * this.im;
	}

	arg() {
		/*
		Computes the angle (argument), as a real number measured in radians
		0 <= arg(z) < 2 * pi
		*/
		return (Math.atan2(this.im, this.re) + 2 * Math.PI) % (2 * Math.PI);
	}

	unit() {
		/* Computes a unit modulus complex number in the direction of this complex number */
		return this.scale(1 / this.norm());
	}

	scale(k) {
		/* Scales each component by the real constant k */
		return new Complex(this.re * k, this.im * k);
	}

	add(z) {
		/* Computes the sum of this complex number and z */
		return new Complex(this.re + z.re, this.im + z.im);
	}

	sub(z) {
		/* Computes the difference of this complex number and z */
		return new Complex(this.re - z.re, this.im - z.im);	
	}

	mult(z) {
		/* Computes the product of this complex number and z */
		return new Complex(this.re * z.re - this.im * z.im, this.re * z.im + this.im * z.re);
	}

	eMult(z) {
		/* elementwise multiplication of this complex number and z */
		return new Complex(this.re * z.re, this.im * z.im);
	}

	inv() {
		/* Computes the reciprocal (inverse) */
		return this.conj().scale(1 / this.normSq());
	}

	div(z) {
		/* Computes the quotient of this complex number and z */
		return this.mult(z.inv());
	}

	perp() {
		/* Computes an orthogonal complex number of the same magnitude */
		return new Complex(-this.im, this.re);
	}

	sqrt() {
		/* Computes the principal branch of the square root */
		const normSqrt = Math.sqrt(this.norm());
		const halfArg = 0.5 * this.arg();
		return new Complex(normSqrt * Math.cos(halfArg), normSqrt * Math.sin(halfArg));
	}

	square() {
		/* Computes the square */
		return new Complex(this.re * this.re - this.im * this.im, 2 * this.re * this.im);
	}

	exp() {
		/* Computes the exponential function of this complex number */
		const mag = Math.exp(this.re);
		return new Complex(mag * Math.cos(this.im), mag * Math.sin(this.im));
	}

	ln() {
		/* Computes the principal branch of the natural log */
		return new Complex(Math.log(this.norm()), this.arg());
	}

	acos() {
		/* Computes the principal branch of the inverse cosine */
		return this.add(this.square().sub(new Complex(1, 0)).sqrt()).ln().div(complex(0, 1));
	}

	asin() {
		/** computes the principle branch of the inverse sine */
		const i = complex(0, 1);
		return this.mult(i).add(complex(1, 0).sub(this.square()).sqrt()).ln().div(i);
	}

	atan() {
		const i = complex(0, 1);
		return Complex.div(Complex.div(i.sub(this), i.add(this)).ln(), i.scale(2));
	}

	rotate(angle) {
		/* Computes this complex number rotated by angle radians */
		return this.mult((new Complex(0, angle)).exp());
	}

	dot(z) {
		/* Computes the Euclidean dot product of the coefficients of this complex number and z */
		return this.re * z.re + this.im * z.im;
	}

	angleTo(z) {
		/* Computes the angle between this complex number and z */
		/*
		acos u*v/uv = uvcos(t)
		*/
		return Math.acos(this.dot(z) / (this.norm() * z.norm()));
	}

	toString() {
		/* Returns the string representation of the complex number as an ordered pair (re(z), im(z)) */
		return `(${this.re},${this.im})`;
	}

	toLatex() {
		/* Returns latex representation of the complex number */
		const rx = roundTo(this.re, 3), ry = roundTo(this.im, 3);
		return `\\left(${rx},${ry}\\right)`
	}

	equals(z) {
		/* Returns true iff z equals this complex number, exactly */
		return (this.re == z.re && this.im == z.im);
	}

	equalsEps(z) {
		/*
		Returns true iff z equals this complex number, within numerical tolerance EPSILON
		For floating point rounding purposes
		*/
		return (Math.abs(this.re - z.re) < EPSILON && Math.abs(this.im - z.im) < EPSILON);
	}

	mobius(a, b, c, d) {
		/*
		Apply the Möbius transformation (az+b)/(cz+d)
		*/
		if (Complex.infinite(this)) {
			if (c !== 0) {
				return Complex.div(a, c);
			} else {
				return new Complex(Infinity, Infinity);
			}
		} else {
			const denominator = Complex.add(Complex.mult(c, this), d);
			if (denominator === 0) {
				return new Complex(Infinity, Infinity);
			} else {
				const numerator = Complex.add(Complex.mult(a, this), b);
				return Complex.div(numerator, denominator);
			}
		}
	}

	sin() {
		/*
		Calculate the sine of the complex number
		*/
		const i = complex(0, 1);
		const rotated = this.mult(i);
		return rotated.exp().sub( rotated.scale(-1).exp() ).div(i.scale(2));
	}

	cos() {
		/*
		Calculate the cosine of the complex number
		*/
		const i = complex(0, 1);
		const rotated = this.mult(i);
		return rotated.exp().add( rotated.scale(-1).exp() ).scale(0.5);
	}

	tan() {
		/*
		Calculate the tangent of the complex number
		*/
		return this.sin().div(this.cos())
	}

	sinh() {
		/*
		Calculate the hyperbolic sine of the complex number
		*/
		return this.exp().sub(this.scale(-1).exp()).scale(0.5);
	}

	cosh() {
		/*
		Calculate the hyperbolic cosine of the complex number
		*/
		return this.exp().add(this.scale(-1).exp()).scale(0.5);
	}

	tanh() {
		/*
		Calculate the hyperbolic tangent of the complex number
		*/
		return this.sinh().div(this.cosh())
	}

	pow(z) {
		/*
		Calculate the zth power of the complex number
		*/
		if (this.equalsEps(complex(0, 0))) {
			return (z.equalsEps(complex(1, 0))) ? complex(1, 0) : complex(0, 0);
		}

		const subAng = Math.atan2(this.im, this.re);
		const normSq = this.normSq();		
		const ang = 0.5 * z.im * Math.log(normSq) + z.re * subAng;
		const norm = Math.exp(-z.im * subAng) * Math.pow(normSq, 0.5 * z.re);
		return complex(norm * Math.cos(ang), norm * Math.sin(ang));
	}

	clamp(min, max) {
		/**
		 * clamp the complex number's norm between two values
		 */
		const norm = z.norm();
		if (!(min <= norm && norm <= max)) {
			return this.scale(clamp(norm, min, max) / norm);
		}
		return this;
	}

	frac() {
		/**
		 * return the fractional part of each of the components of the complex number
		 */
		return complex(
			frac(this.re),
			frac(this.im),
		);
	}

	/* ---------- Static functions -------------------- */

	static norm(z) {
		return z.norm();
	}

	static normSq(z) {
		return z.normSq();
	}

	static inv(z) {
		return z.conj().scale(1 / z.normSq());
	}

	static scale(z, k) {
		return z.scale(k);
	}

	static add(z1, z2) {
		return z1.add(z2);
	}

	static sub(z1, z2) {
		return z1.sub(z2);
	}

	static mult(z1, z2) {
		return z1.mult(z2);
	}

	static div(z1, z2) {
		return z1.div(z2);
	}

	static pow(z1, z2) {
		return z1.pow(z2);
	}

	static exp(z) {
		return z.exp();
	}

	static ln(z) {
		return z.ln();
	}

	static sqrt(z) {
		return z.sqrt();
	}

	static sin(z) {
		/*
		Return the sine of z
		*/
		return z.sin();
	}

	static cos(z) {
		/*
		Return the cosine of z
		*/
		return z.cos();
	}

	static tan(z) {
		/*
		Return the tangent of z		
		*/
		return z.tan();
	}

	static asin(z) {
		return z.asin();
	} 

	static acos(z) {
		return z.acos();
	}

	static atan(z) {
		return z.atan();
	}

	static sinh(z) {
		/*
		Return the hyperbolic sine of z
		*/
		return z.sinh();
	}

	static cosh(z) {
		/*
		Return the hyperbolic cosine of z
		*/
		return z.cosh();
	}

	static tanh(z) {
		/*
		Return the hyperbolic tangent of z
		*/
		return z.tanh();
	}

	static infinite(z) {
		/**
		 * return whether one or both of the components of z is infinite, or if the norm is infinite
		 */
		const norm = z.norm();
		return (
			z.re === Infinity || z.re === -Infinity 
			|| z.im === Infinity || z.im === -Infinity
			|| norm === Infinity
		);
	}

	static nan(z) {
		/**
		 * return whether one or both components of z is NaN 
		*/
		return isNaN(z.re) || isNaN(z.im);
	}

	static gamma(z) {
		if (z.re < 0.5) {
			// gamma(1-z)gamma(z) = pi / sin(pi*z)
			return Complex.div( complex(Math.PI, 0.0), Complex.mult(z.scale(Math.PI).sin(), Complex.gamma(complex(1 - z.re, -z.im))) );
		} else {
			z = Complex.sub(z, complex(1, 0)); // account for stupid shift by 1
			let x = complex(pValues[0], 0);
			for (let i=1; i<pValues.length; i++) {
				x = Complex.add(x, Complex.div( complex(pValues[i], 0.0), Complex.add(z, complex(i, 0)) ));
			}
			const t = Complex.add(z, complex(7.5, 0)); // g=7, g+0.5
			return Complex.pow(t, z.add(complex(0.5, 0))).mult(t.scale(-1).exp()).mult(x).scale(Math.sqrt(2 * Math.PI));
		}
	}

	static beta(z1, z2) {
		return Complex.div(Complex.mult(Complex.gamma(z1), Complex.gamma(z2)), Complex.gamma(z1.add(z2)));
	}

	static max(z1, z2) {
		return (z1.normSq() < z2.normSq()) ? z2 : z1;
	}

	static min(z1, z2) {
		return (z1.normSq() > z2.normSq()) ? z2 : z1;
	}

	static clamp(z, min, max) {
		return z.clamp(min, max);
	}

	static frac() {
		return z.frac();
	}

	/* --------------- In-place operations --------------------- */

	iadd(z) {
		this.re += z.re;
		this.im += z.im;
	}

	isub(z) {
		this.re -= z.re;
		this.im -= z.im;
	}


}

function complex(real, imaginary) {
	/* instantiate a Complex without new keyword */
	return new Complex(real, imaginary);
}


function integrateOverParameter(f, a=0, b=1, samples=100) {
	/*
	return the line integral of the Complex-valued function f
	parameterized by the real parameter t on [a, b],
	using the specified number of samples
	*/
	const step = (b - a) / samples;
	let t = a + step / 2; // midpoint approximation
	let result = complex(0, 0);
	for (let i=0; i<samples; i++) {
		result.iadd(f(t).scale(step));
		t += step;
	}
	return result;
}

function integrateOverComplexVariable(f, z0, z1, samples=100) {
	/*
	return the line integral of the Complex-valued function f
	on the line between z0 and z1, using the specified number of intervals
	*/
	const step = Complex.sub(z1, z0).scale(1 / samples);
	const stepNorm = step.norm();
	let z = Complex.add(z0, step.scale(0.5));
	let result = complex(0, 0);
	for (let i=0; i<samples; i++) {
		result.iadd(f(z).scale(stepNorm));
		z.iadd(step);
	}
	return result;
}



class Vector {

	/*
	Represents a vector in C^n
	*/

	constructor(values) {
		this.values = values;
		this.dimension = values.length;
	}

	get(index) {
		return this.values[index];
	}

	scale(z) {
		/*
		Scale each component by the complex number z
		*/
		const result = [];
		for (let i=0; i<this.dimension; i++) {
			result.push(Complex.mult(this.get(i), z));
		}
		return vector(result);
	}

	realScale(k) {
		/*
		Scale each component by the real number k
		*/
		const result = [];
		for (let i=0; i<this.dimension; i++) {
			result.push(this.get(i).scale(k));
		}
		return vector(result);
	}

	add(Z) {
		const result = [];
		for (let i=0; i<this.dimension; i++) {
			result.push(Complex.add(this.get(i), Z.get(i)));
		}
		return vector(result);
	}
	
	dot(Z) {
		/*
		Compute (Hermitian) dot product
		Equivalent to standard dot product in the case of real vectors
		*/
		const result = complex(0, 0);
		for (let i=0; i<this.dimension; i++) {
			result.iadd(Complex.mult(this.get(i), Z.get(i).conj()));
		}
		return result;
	}

	norm() {
		/*
		Compute norm of the vector
		*/
		return Math.sqrt(this.dot(this).re);
	}

	angleTo(Z) {
		/*
		Compute the angle between the vector and vector Z
		*/
		return Math.acos( this.dot(Z).re / (this.norm() * Z.norm()) );
	}

	static lerp(Z1, Z2, t) {
		/*
		Linearly interpolate
		*/
		return Z1.realScale(1 - t).add(Z2.scale(t));
	}

}


function vector(...values) {
	/*
	Utility function for easy instantiation of Vector objects.
	Allows components to be passed individually or in an array
	*/
	if (values[0] instanceof Complex) {
		return new Vector(values);
	} else {
		return new Vector(values[0]);
	}
}


module.exports = {
	Complex, complex,
	integrateOverComplexVariable, integrateOverParameter,
	Vector, vector,
};
},{}]},{},[2])
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uL1VzZXJzL3lhbmtlL0FwcERhdGEvUm9hbWluZy9ucG0vbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImFwcC9zY29wZS5qcyIsImRvY3MvZ3VpZGVfZ2VuLmpzIiwibWF0aC9jb21wbGV4LmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0UUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JEQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24oKXtmdW5jdGlvbiByKGUsbix0KXtmdW5jdGlvbiBvKGksZil7aWYoIW5baV0pe2lmKCFlW2ldKXt2YXIgYz1cImZ1bmN0aW9uXCI9PXR5cGVvZiByZXF1aXJlJiZyZXF1aXJlO2lmKCFmJiZjKXJldHVybiBjKGksITApO2lmKHUpcmV0dXJuIHUoaSwhMCk7dmFyIGE9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitpK1wiJ1wiKTt0aHJvdyBhLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsYX12YXIgcD1uW2ldPXtleHBvcnRzOnt9fTtlW2ldWzBdLmNhbGwocC5leHBvcnRzLGZ1bmN0aW9uKHIpe3ZhciBuPWVbaV1bMV1bcl07cmV0dXJuIG8obnx8cil9LHAscC5leHBvcnRzLHIsZSxuLHQpfXJldHVybiBuW2ldLmV4cG9ydHN9Zm9yKHZhciB1PVwiZnVuY3Rpb25cIj09dHlwZW9mIHJlcXVpcmUmJnJlcXVpcmUsaT0wO2k8dC5sZW5ndGg7aSsrKW8odFtpXSk7cmV0dXJuIG99cmV0dXJuIHJ9KSgpIiwiXHJcbmNvbnN0IHsgY29tcGxleCwgQ29tcGxleCB9ID0gcmVxdWlyZShcIi4uL21hdGgvY29tcGxleC5qc1wiKTtcclxuXHJcblxyXG5jb25zdCBzY29wZSA9IHtcclxuICAgIGJ1aWx0aW46IHtcclxuICAgICAgICBcInpcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwielwiLFxyXG4gICAgICAgICAgICBpc1BhcmFtZXRlcjogdHJ1ZSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIC8vIFwidFwiOiB7XHJcbiAgICAgICAgLy8gICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgIC8vICAgICBzaGFkZXJBbGlhczogXCJ0XCIsXHJcbiAgICAgICAgLy8gICAgIGlzUGFyYW1ldGVyOiB0cnVlLFxyXG4gICAgICAgIC8vIH0sXHJcblxyXG4gICAgICAgIFwibm9ybVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcIm5vcm1DXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibm9ybVNxXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibm9ybVNxQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImFyZ1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImFyZ0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJpbnZcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpbnZDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwiZXhwXCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImV4cENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJsblwiOiB7IFxyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJsbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJzcXJ0XCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInNxcnRDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwic2luXCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInNpbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJjb3NcIjogeyBcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY29zQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInRhblwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInRhbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhc2luXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYXNpbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhY29zXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYWNvc0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhdGFuXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYXRhbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJzaW5oXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwic2luaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJjb3NoXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY29zaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJ0YW5oXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwidGFuaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhdGFuaFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImF0YW5oQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInJlXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwicmVDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwiaW1cIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpbUNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJHYW1tYVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcIkdhbW1hQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImJldGFcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJiZXRhQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibWluXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibWluQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibWF4XCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibWF4Q1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibGVycFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImxlcnBDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwid1wiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSwgXCJ0XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMiB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImNvbmpcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJjb25qQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImNsYW1wXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY2xhbXBDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwid1wiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSwgXCJ0XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMiB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImZyYWNcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJmcmFjQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImludmVyc2VTQ1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImludmVyc2VTQ0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSwgXCJwXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInNjXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwic2NDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwicFwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJwbGFuZVRvUFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInBsYW5lVG9QQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcInBcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwicFRvUGxhbmVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJwVG9QbGFuZUNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSwgXCJwXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMSB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInNxdWVlemVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJzcXVlZXplQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcImNvdmVyYWdlXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMSB9LCBcImxlbmd0aFwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDIgfSB9LFxyXG4gICAgICAgIH0sXHJcblxyXG4gICAgICAgIFwiaVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpXCIsXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInBpXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogZmFsc2UsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInBpQ1wiLFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJ0YXVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwidHBpQ1wiLFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJlXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogZmFsc2UsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImVcIixcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwicmVhbEJvdW5kc1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJ4Qm91bmRzXCIsXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImltYWdCb3VuZHNcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwieUJvdW5kc1wiLFxyXG4gICAgICAgIH0sXHJcblxyXG4gICAgICAgIC8qIGRhdGF0eXBlcyAoZm9yIHR5cGUgYW5ub3RhdGlvbikgKi9cclxuICAgICAgICBcImNvbXBsZXhcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgaXNUeXBlOiB0cnVlLCB9LFxyXG4gICAgICAgIFwiYXJyYXlcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgaXNUeXBlOiB0cnVlLCB9LFxyXG4gICAgICAgIFwibWF0cml4XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIGlzVHlwZTogdHJ1ZSwgfSxcclxuICAgICAgICBcImZ1bmN0aW9uXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIGlzVHlwZTogdHJ1ZSwgfSxcclxuICAgIH0sXHJcblxyXG4gICAgdXNlckdsb2JhbDoge1xyXG4gICAgfSxcclxufVxyXG5cclxuXHJcbmNvbnN0IGRlZmF1bHRWYWx1ZVNjb3BlID0ge1xyXG4gICAgXCJ6XCI6IGNvbXBsZXgoMSwgMCksXHJcbiAgICBcIm5vcm1cIjogKHopID0+IGNvbXBsZXgoei5ub3JtKCksIDApLFxyXG4gICAgXCJub3JtU3FcIjogKHopID0+IGNvbXBsZXgoei5ub3JtU3EoKSwgMCksXHJcbiAgICBcImFyZ1wiOiAoeikgPT4gY29tcGxleCh6LmFyZygpLCAwKSxcclxuICAgIFwiaW52XCI6IENvbXBsZXguaW52LFxyXG4gICAgXCJleHBcIjogQ29tcGxleC5leHAsXHJcbiAgICBcImxuXCI6IENvbXBsZXgubG4sXHJcbiAgICBcInNxcnRcIjogQ29tcGxleC5zcXJ0LFxyXG4gICAgXCJzaW5cIjogQ29tcGxleC5zaW4sXHJcbiAgICBcImNvc1wiOiBDb21wbGV4LmNvcyxcclxuICAgIFwidGFuXCI6IENvbXBsZXgudGFuLFxyXG4gICAgXCJhY29zXCI6IENvbXBsZXguYWNvcyxcclxuICAgIFwiYXNpblwiOiBDb21wbGV4LmFzaW4sXHJcbiAgICBcImF0YW5cIjogQ29tcGxleC5hdGFuLFxyXG4gICAgXCJzaW5oXCI6IENvbXBsZXguc2luaCxcclxuICAgIFwiY29zaFwiOiBDb21wbGV4LmNvc2gsXHJcbiAgICBcInRhbmhcIjogQ29tcGxleC50YW5oLFxyXG4gICAgXCJhdGFuaFwiOiAoeikgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInJlXCI6ICh6KSA9PiBjb21wbGV4KHoucmUsIDApLFxyXG4gICAgXCJpbVwiOiAoeikgPT4gY29tcGxleCh6LmltLCAwKSxcclxuICAgIFwiR2FtbWFcIjogQ29tcGxleC5nYW1tYSxcclxuICAgIFwiYmV0YVwiOiBDb21wbGV4LmJldGEsXHJcbiAgICBcIm1pblwiOiBDb21wbGV4Lm1pbixcclxuICAgIFwibWF4XCI6IENvbXBsZXgubWF4LFxyXG4gICAgXCJsZXJwXCI6ICh6MSwgejIsIHQpID0+IENvbXBsZXgubXVsdChjb21wbGV4KDEsIDApLnN1Yih0KSwgejEpLmFkZCh6Mi5tdWx0KHQpKSxcclxuICAgIFwiY29ualwiOiAoeikgPT4gei5jb25qKCksXHJcbiAgICBcImNsYW1wXCI6ICh6LCBtaW4sIG1heCkgPT4gQ29tcGxleC5jbGFtcCh6LCBtaW4ubm9ybSgpLCBtYXgubm9ybSgpKSxcclxuICAgIFwiZnJhY1wiOiBDb21wbGV4LmZyYWMsXHJcbiAgICBcImludmVyc2VTQ1wiOiAoeiwgcCkgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInNjXCI6ICh6LCBwKSA9PiBjb21wbGV4KDEsIDApLCAvLyBub3QgaW1wbGVtZW50ZWQgZm9yIHA1IG1vZGVcclxuICAgIFwicGxhbmVUb1BcIjogKHosIHApID0+IGNvbXBsZXgoMSwgMCksIC8vIG5vdCBpbXBsZW1lbnRlZCBmb3IgcDUgbW9kZVxyXG4gICAgXCJwVG9QbGFuZVwiOiAoeiwgcCkgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInNxdWVlemVcIjogKHosIGNvdmVyYWdlLCBsZW5ndGgpID0+IGNvbXBsZXgoMSwgMCksIC8vIG5vdCBpbXBsZW1lbnRlZCBmb3IgcDUgbW9kZVxyXG4gIFxyXG4gICAgXCJpXCI6IGNvbXBsZXgoMCwgMSksXHJcbiAgICBcInBpXCI6IGNvbXBsZXgoTWF0aC5QSSwgMCksXHJcbiAgICBcInRhdVwiOiBjb21wbGV4KDIgKiBNYXRoLlBJLCAwKSxcclxuICAgIFwiZVwiOiBjb21wbGV4KE1hdGguRSwgMCksXHJcbiAgICBcInJlYWxCb3VuZHNcIjogY29tcGxleCgwLCAwKSxcclxuICAgIFwiaW1hZ0JvdW5kc1wiOiBjb21wbGV4KDAsIDApLFxyXG59O1xyXG5cclxubGV0IHZhbHVlU2NvcGUgPSB7fTtcclxuXHJcbm1vZHVsZS5leHBvcnRzID0ge1xyXG4gICAgc2NvcGUsIGRlZmF1bHRWYWx1ZVNjb3BlLCB2YWx1ZVNjb3BlLFxyXG59OyIsIlxyXG5cclxuXHJcbmNvbnN0IHNjb3BlcyA9IHJlcXVpcmUoXCIuLi9hcHAvc2NvcGUuanNcIik7XHJcbmNvbnN0IHNjb3BlID0gc2NvcGVzLnNjb3BlLmJ1aWx0aW47XHJcblxyXG5jb25zdCBkZXNjcmlwdGlvbnMgPSB7XHJcbiAgICBcIkdhbW1hXCI6IFwiVGhlIEdhbW1hIGZ1bmN0aW9uIFxcXFwoIFxcXFxHYW1tYSh6KSBcXFxcKVwiLFxyXG59O1xyXG5cclxuXHJcbmNvbnN0IGdlbmVyYXRlRGVzY3JpcHRpb25zID0gKCkgPT4ge1xyXG4gICAgY29uc3Qga2V5cyA9IE9iamVjdC5rZXlzKHNjb3BlKS5zb3J0KCk7XHJcbiAgICBsZXQgcmVzdWx0ID0gXCJcIjtcclxuXHJcbiAgICBmb3IgKGNvbnN0IGtleSBvZiBrZXlzKSB7XHJcbiAgICAgICAgY29uc3QgZGVzY3JpcHRpb24gPSBkZXNjcmlwdGlvbnNba2V5XSA/PyBcIk5vIGRlc2NyaXB0aW9uIGdpdmVuXCI7XHJcbiAgICAgICAgbGV0IGFyZ3VtZW50cyA9IFwiXCI7XHJcbiAgICAgICAgaWYgKHNjb3BlW2tleV0uaXNGdW5jdGlvbikge1xyXG4gICAgICAgICAgICBhcmd1bWVudHMgPSBcIkFyZ3VtZW50czogXCI7XHJcbiAgICAgICAgICAgIGNvbnN0IGxvY2FscyA9IE9iamVjdC5rZXlzKHNjb3BlW2tleV0ubG9jYWxzKS5zb3J0KGxvY2FsID0+IGxvY2FsLmluZGV4KTtcclxuICAgICAgICAgICAgaWYgKGxvY2Fscy5sZW5ndGggPT09IDApIGFyZ3VtZW50cyArPSBcIm5vbmVcIjtcclxuICAgICAgICAgICAgZm9yIChsZXQgaT0wOyBpPGxvY2Fscy5sZW5ndGg7IGkrKykge1xyXG4gICAgICAgICAgICAgICAgY29uc3QgbG9jYWwgPSBsb2NhbHNbaV07XHJcbiAgICAgICAgICAgICAgICBhcmd1bWVudHMgKz0gbG9jYWwgKyBgICgke3Njb3BlW2tleV0ubG9jYWxzW2xvY2FsXS50eXBlfSlgO1xyXG4gICAgICAgICAgICAgICAgaWYgKGkgIT09IGxvY2Fscy5sZW5ndGggLSAxKSBhcmd1bWVudHMgKz0gXCIsIFwiO1xyXG4gICAgICAgICAgICB9XHJcbiAgICAgICAgICAgIGFyZ3VtZW50cyArPSBcIjxicj5cIjtcclxuICAgICAgICB9XHJcbiAgICAgICAgcmVzdWx0ICs9IGA8ZGl2IGNsYXNzPVwiZGVzY3JpcHRpb24tZW50cnlcIj48c3BhbiBzdHlsZT1cImZvbnQtd2VpZ2h0OiBib2xkXCI+JHtrZXl9PC9zcGFuPjxicj4ke2FyZ3VtZW50c31EZXNjcmlwdGlvbjoke2Rlc2NyaXB0aW9ufTwvZGl2PmA7XHJcbiAgICB9XHJcblxyXG4gICAgcmV0dXJuIHJlc3VsdDtcclxufTtcclxuXHJcbmNvbnN0IG91dHB1dERlc2NyaXB0aW9ucyA9ICh0YXJnZXRJZCkgPT4ge1xyXG4gICAgY29uc3QgdGFyZ2V0RWxlbWVudCA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoYCMke3RhcmdldElkfWApO1xyXG4gICAgdGFyZ2V0RWxlbWVudC5pbm5lckhUTUwgPSBnZW5lcmF0ZURlc2NyaXB0aW9ucygpO1xyXG59XHJcblxyXG5jb25zdCB0b2dnbGVEZXNjcmlwdGlvbnMgPSAoKSA9PiB7XHJcbiAgICBjb25zdCBkZXNjRGl2ID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihcIiNkZXNjcmlwdGlvbi1jb250YWluZXJcIik7XHJcbiAgICBjb25zdCBjb2xsYXBzZUJ0biA9IGRvY3VtZW50LnF1ZXJ5U2VsZWN0b3IoXCIjY29sbGFwc2UtYnRuXCIpO1xyXG4gICAgaWYgKGRlc2NEaXYuc3R5bGUuZGlzcGxheSA9PT0gXCJub25lXCIpIHtcclxuICAgICAgICBkZXNjRGl2LnN0eWxlLmRpc3BsYXkgPSBcImJsb2NrXCI7XHJcbiAgICAgICAgY29sbGFwc2VCdG4uaW5uZXJUZXh0ID0gXCJDb2xsYXBzZVwiO1xyXG4gICAgfSBlbHNlIHtcclxuICAgICAgICBkZXNjRGl2LnN0eWxlLmRpc3BsYXkgPSBcIm5vbmVcIjtcclxuICAgICAgICBjb2xsYXBzZUJ0bi5pbm5lclRleHQgPSBcIkV4cGFuZFwiO1xyXG4gICAgfVxyXG59XHJcblxyXG5vdXRwdXREZXNjcmlwdGlvbnMoXCJkZXNjcmlwdGlvbi1jb250YWluZXJcIik7XHJcbndpbmRvdy50b2dnbGVEZXNjcmlwdGlvbnMgPSB0b2dnbGVEZXNjcmlwdGlvbnM7IiwiY29uc3QgRVBTSUxPTiA9IDAuMDAwMDAxO1xyXG5cclxuY29uc3QgcFZhbHVlcyA9IFtcclxuXHQwLjk5OTk5OTk5OTk5OTgwOTkzLFxyXG5cdDY3Ni41MjAzNjgxMjE4ODUxLFxyXG5cdC0xMjU5LjEzOTIxNjcyMjQwMjgsXHJcblx0NzcxLjMyMzQyODc3NzY1MzEzLFxyXG5cdC0xNzYuNjE1MDI5MTYyMTQwNTksXHJcblx0MTIuNTA3MzQzMjc4Njg2OTA1LFxyXG5cdC0wLjEzODU3MTA5NTI2NTcyMDEyLFxyXG5cdDkuOTg0MzY5NTc4MDE5NTcxNmUtNixcclxuXHQxLjUwNTYzMjczNTE0OTMxMTZlLTdcclxuXTtcclxuXHJcblxyXG5jb25zdCBmcmFjID0gKHgpID0+IHtcclxuXHRyZXR1cm4geCAtICh4ID4gMCA/IE1hdGguZmxvb3IoeCkgOiBNYXRoLmNlaWwoeCkpO1xyXG59O1xyXG5cclxuY29uc3QgY2xhbXAgPSAoeCwgbWluLCBtYXgpID0+IHtcclxuXHRyZXR1cm4gTWF0aC5taW4obWF4LCBNYXRoLm1heChtaW4sIHgpKTtcclxufTtcclxuXHJcblxyXG5jbGFzcyBDb21wbGV4IHtcclxuXHJcblx0LypcclxuXHRDbGFzcyBmb3IgcmVwcmVzZW50aW5nIGNvbXBsZXggbnVtYmVycyBvZiB0aGUgZm9ybSBhICsgYmlcclxuXHQqL1xyXG5cclxuXHRjb25zdHJ1Y3RvcihyZWFsLCBpbWFnaW5hcnkpIHtcclxuXHRcdHRoaXMucmUgPSByZWFsO1xyXG5cdFx0dGhpcy5pbSA9IGltYWdpbmFyeTtcclxuXHR9XHJcblxyXG5cdGNvbmooKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgY29tcGxleCBjb25qdWdhdGUgKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleCh0aGlzLnJlLCAtdGhpcy5pbSk7XHJcblx0fVxyXG5cclxuXHRub3JtKCkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIG5vcm0gKG1vZHVsdXMpLCBhcyBhIHJlYWwgbnVtYmVyICovXHJcblx0XHRyZXR1cm4gTWF0aC5zcXJ0KHRoaXMucmUgKiB0aGlzLnJlICsgdGhpcy5pbSAqIHRoaXMuaW0pO1xyXG5cdH1cclxuXHJcblx0bm9ybVNxKCkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIHNxdWFyZSBvZiB0aGUgbm9ybSAobW9kdWx1cyksIGFzIGEgcmVhbCBudW1iZXIgKi9cclxuXHRcdHJldHVybiB0aGlzLnJlICogdGhpcy5yZSArIHRoaXMuaW0gKiB0aGlzLmltO1xyXG5cdH1cclxuXHJcblx0YXJnKCkge1xyXG5cdFx0LypcclxuXHRcdENvbXB1dGVzIHRoZSBhbmdsZSAoYXJndW1lbnQpLCBhcyBhIHJlYWwgbnVtYmVyIG1lYXN1cmVkIGluIHJhZGlhbnNcclxuXHRcdDAgPD0gYXJnKHopIDwgMiAqIHBpXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIChNYXRoLmF0YW4yKHRoaXMuaW0sIHRoaXMucmUpICsgMiAqIE1hdGguUEkpICUgKDIgKiBNYXRoLlBJKTtcclxuXHR9XHJcblxyXG5cdHVuaXQoKSB7XHJcblx0XHQvKiBDb21wdXRlcyBhIHVuaXQgbW9kdWx1cyBjb21wbGV4IG51bWJlciBpbiB0aGUgZGlyZWN0aW9uIG9mIHRoaXMgY29tcGxleCBudW1iZXIgKi9cclxuXHRcdHJldHVybiB0aGlzLnNjYWxlKDEgLyB0aGlzLm5vcm0oKSk7XHJcblx0fVxyXG5cclxuXHRzY2FsZShrKSB7XHJcblx0XHQvKiBTY2FsZXMgZWFjaCBjb21wb25lbnQgYnkgdGhlIHJlYWwgY29uc3RhbnQgayAqL1xyXG5cdFx0cmV0dXJuIG5ldyBDb21wbGV4KHRoaXMucmUgKiBrLCB0aGlzLmltICogayk7XHJcblx0fVxyXG5cclxuXHRhZGQoeikge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIHN1bSBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSArIHoucmUsIHRoaXMuaW0gKyB6LmltKTtcclxuXHR9XHJcblxyXG5cdHN1Yih6KSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgZGlmZmVyZW5jZSBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSAtIHoucmUsIHRoaXMuaW0gLSB6LmltKTtcdFxyXG5cdH1cclxuXHJcblx0bXVsdCh6KSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgcHJvZHVjdCBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSAqIHoucmUgLSB0aGlzLmltICogei5pbSwgdGhpcy5yZSAqIHouaW0gKyB0aGlzLmltICogei5yZSk7XHJcblx0fVxyXG5cclxuXHRlTXVsdCh6KSB7XHJcblx0XHQvKiBlbGVtZW50d2lzZSBtdWx0aXBsaWNhdGlvbiBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSAqIHoucmUsIHRoaXMuaW0gKiB6LmltKTtcclxuXHR9XHJcblxyXG5cdGludigpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSByZWNpcHJvY2FsIChpbnZlcnNlKSAqL1xyXG5cdFx0cmV0dXJuIHRoaXMuY29uaigpLnNjYWxlKDEgLyB0aGlzLm5vcm1TcSgpKTtcclxuXHR9XHJcblxyXG5cdGRpdih6KSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgcXVvdGllbnQgb2YgdGhpcyBjb21wbGV4IG51bWJlciBhbmQgeiAqL1xyXG5cdFx0cmV0dXJuIHRoaXMubXVsdCh6LmludigpKTtcclxuXHR9XHJcblxyXG5cdHBlcnAoKSB7XHJcblx0XHQvKiBDb21wdXRlcyBhbiBvcnRob2dvbmFsIGNvbXBsZXggbnVtYmVyIG9mIHRoZSBzYW1lIG1hZ25pdHVkZSAqL1xyXG5cdFx0cmV0dXJuIG5ldyBDb21wbGV4KC10aGlzLmltLCB0aGlzLnJlKTtcclxuXHR9XHJcblxyXG5cdHNxcnQoKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgcHJpbmNpcGFsIGJyYW5jaCBvZiB0aGUgc3F1YXJlIHJvb3QgKi9cclxuXHRcdGNvbnN0IG5vcm1TcXJ0ID0gTWF0aC5zcXJ0KHRoaXMubm9ybSgpKTtcclxuXHRcdGNvbnN0IGhhbGZBcmcgPSAwLjUgKiB0aGlzLmFyZygpO1xyXG5cdFx0cmV0dXJuIG5ldyBDb21wbGV4KG5vcm1TcXJ0ICogTWF0aC5jb3MoaGFsZkFyZyksIG5vcm1TcXJ0ICogTWF0aC5zaW4oaGFsZkFyZykpO1xyXG5cdH1cclxuXHJcblx0c3F1YXJlKCkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIHNxdWFyZSAqL1xyXG5cdFx0cmV0dXJuIG5ldyBDb21wbGV4KHRoaXMucmUgKiB0aGlzLnJlIC0gdGhpcy5pbSAqIHRoaXMuaW0sIDIgKiB0aGlzLnJlICogdGhpcy5pbSk7XHJcblx0fVxyXG5cclxuXHRleHAoKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgZXhwb25lbnRpYWwgZnVuY3Rpb24gb2YgdGhpcyBjb21wbGV4IG51bWJlciAqL1xyXG5cdFx0Y29uc3QgbWFnID0gTWF0aC5leHAodGhpcy5yZSk7XHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgobWFnICogTWF0aC5jb3ModGhpcy5pbSksIG1hZyAqIE1hdGguc2luKHRoaXMuaW0pKTtcclxuXHR9XHJcblxyXG5cdGxuKCkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIHByaW5jaXBhbCBicmFuY2ggb2YgdGhlIG5hdHVyYWwgbG9nICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgoTWF0aC5sb2codGhpcy5ub3JtKCkpLCB0aGlzLmFyZygpKTtcclxuXHR9XHJcblxyXG5cdGFjb3MoKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgcHJpbmNpcGFsIGJyYW5jaCBvZiB0aGUgaW52ZXJzZSBjb3NpbmUgKi9cclxuXHRcdHJldHVybiB0aGlzLmFkZCh0aGlzLnNxdWFyZSgpLnN1YihuZXcgQ29tcGxleCgxLCAwKSkuc3FydCgpKS5sbigpLmRpdihjb21wbGV4KDAsIDEpKTtcclxuXHR9XHJcblxyXG5cdGFzaW4oKSB7XHJcblx0XHQvKiogY29tcHV0ZXMgdGhlIHByaW5jaXBsZSBicmFuY2ggb2YgdGhlIGludmVyc2Ugc2luZSAqL1xyXG5cdFx0Y29uc3QgaSA9IGNvbXBsZXgoMCwgMSk7XHJcblx0XHRyZXR1cm4gdGhpcy5tdWx0KGkpLmFkZChjb21wbGV4KDEsIDApLnN1Yih0aGlzLnNxdWFyZSgpKS5zcXJ0KCkpLmxuKCkuZGl2KGkpO1xyXG5cdH1cclxuXHJcblx0YXRhbigpIHtcclxuXHRcdGNvbnN0IGkgPSBjb21wbGV4KDAsIDEpO1xyXG5cdFx0cmV0dXJuIENvbXBsZXguZGl2KENvbXBsZXguZGl2KGkuc3ViKHRoaXMpLCBpLmFkZCh0aGlzKSkubG4oKSwgaS5zY2FsZSgyKSk7XHJcblx0fVxyXG5cclxuXHRyb3RhdGUoYW5nbGUpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoaXMgY29tcGxleCBudW1iZXIgcm90YXRlZCBieSBhbmdsZSByYWRpYW5zICovXHJcblx0XHRyZXR1cm4gdGhpcy5tdWx0KChuZXcgQ29tcGxleCgwLCBhbmdsZSkpLmV4cCgpKTtcclxuXHR9XHJcblxyXG5cdGRvdCh6KSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgRXVjbGlkZWFuIGRvdCBwcm9kdWN0IG9mIHRoZSBjb2VmZmljaWVudHMgb2YgdGhpcyBjb21wbGV4IG51bWJlciBhbmQgeiAqL1xyXG5cdFx0cmV0dXJuIHRoaXMucmUgKiB6LnJlICsgdGhpcy5pbSAqIHouaW07XHJcblx0fVxyXG5cclxuXHRhbmdsZVRvKHopIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBhbmdsZSBiZXR3ZWVuIHRoaXMgY29tcGxleCBudW1iZXIgYW5kIHogKi9cclxuXHRcdC8qXHJcblx0XHRhY29zIHUqdi91diA9IHV2Y29zKHQpXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIE1hdGguYWNvcyh0aGlzLmRvdCh6KSAvICh0aGlzLm5vcm0oKSAqIHoubm9ybSgpKSk7XHJcblx0fVxyXG5cclxuXHR0b1N0cmluZygpIHtcclxuXHRcdC8qIFJldHVybnMgdGhlIHN0cmluZyByZXByZXNlbnRhdGlvbiBvZiB0aGUgY29tcGxleCBudW1iZXIgYXMgYW4gb3JkZXJlZCBwYWlyIChyZSh6KSwgaW0oeikpICovXHJcblx0XHRyZXR1cm4gYCgke3RoaXMucmV9LCR7dGhpcy5pbX0pYDtcclxuXHR9XHJcblxyXG5cdHRvTGF0ZXgoKSB7XHJcblx0XHQvKiBSZXR1cm5zIGxhdGV4IHJlcHJlc2VudGF0aW9uIG9mIHRoZSBjb21wbGV4IG51bWJlciAqL1xyXG5cdFx0Y29uc3QgcnggPSByb3VuZFRvKHRoaXMucmUsIDMpLCByeSA9IHJvdW5kVG8odGhpcy5pbSwgMyk7XHJcblx0XHRyZXR1cm4gYFxcXFxsZWZ0KCR7cnh9LCR7cnl9XFxcXHJpZ2h0KWBcclxuXHR9XHJcblxyXG5cdGVxdWFscyh6KSB7XHJcblx0XHQvKiBSZXR1cm5zIHRydWUgaWZmIHogZXF1YWxzIHRoaXMgY29tcGxleCBudW1iZXIsIGV4YWN0bHkgKi9cclxuXHRcdHJldHVybiAodGhpcy5yZSA9PSB6LnJlICYmIHRoaXMuaW0gPT0gei5pbSk7XHJcblx0fVxyXG5cclxuXHRlcXVhbHNFcHMoeikge1xyXG5cdFx0LypcclxuXHRcdFJldHVybnMgdHJ1ZSBpZmYgeiBlcXVhbHMgdGhpcyBjb21wbGV4IG51bWJlciwgd2l0aGluIG51bWVyaWNhbCB0b2xlcmFuY2UgRVBTSUxPTlxyXG5cdFx0Rm9yIGZsb2F0aW5nIHBvaW50IHJvdW5kaW5nIHB1cnBvc2VzXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIChNYXRoLmFicyh0aGlzLnJlIC0gei5yZSkgPCBFUFNJTE9OICYmIE1hdGguYWJzKHRoaXMuaW0gLSB6LmltKSA8IEVQU0lMT04pO1xyXG5cdH1cclxuXHJcblx0bW9iaXVzKGEsIGIsIGMsIGQpIHtcclxuXHRcdC8qXHJcblx0XHRBcHBseSB0aGUgTcO2Yml1cyB0cmFuc2Zvcm1hdGlvbiAoYXorYikvKGN6K2QpXHJcblx0XHQqL1xyXG5cdFx0aWYgKENvbXBsZXguaW5maW5pdGUodGhpcykpIHtcclxuXHRcdFx0aWYgKGMgIT09IDApIHtcclxuXHRcdFx0XHRyZXR1cm4gQ29tcGxleC5kaXYoYSwgYyk7XHJcblx0XHRcdH0gZWxzZSB7XHJcblx0XHRcdFx0cmV0dXJuIG5ldyBDb21wbGV4KEluZmluaXR5LCBJbmZpbml0eSk7XHJcblx0XHRcdH1cclxuXHRcdH0gZWxzZSB7XHJcblx0XHRcdGNvbnN0IGRlbm9taW5hdG9yID0gQ29tcGxleC5hZGQoQ29tcGxleC5tdWx0KGMsIHRoaXMpLCBkKTtcclxuXHRcdFx0aWYgKGRlbm9taW5hdG9yID09PSAwKSB7XHJcblx0XHRcdFx0cmV0dXJuIG5ldyBDb21wbGV4KEluZmluaXR5LCBJbmZpbml0eSk7XHJcblx0XHRcdH0gZWxzZSB7XHJcblx0XHRcdFx0Y29uc3QgbnVtZXJhdG9yID0gQ29tcGxleC5hZGQoQ29tcGxleC5tdWx0KGEsIHRoaXMpLCBiKTtcclxuXHRcdFx0XHRyZXR1cm4gQ29tcGxleC5kaXYobnVtZXJhdG9yLCBkZW5vbWluYXRvcik7XHJcblx0XHRcdH1cclxuXHRcdH1cclxuXHR9XHJcblxyXG5cdHNpbigpIHtcclxuXHRcdC8qXHJcblx0XHRDYWxjdWxhdGUgdGhlIHNpbmUgb2YgdGhlIGNvbXBsZXggbnVtYmVyXHJcblx0XHQqL1xyXG5cdFx0Y29uc3QgaSA9IGNvbXBsZXgoMCwgMSk7XHJcblx0XHRjb25zdCByb3RhdGVkID0gdGhpcy5tdWx0KGkpO1xyXG5cdFx0cmV0dXJuIHJvdGF0ZWQuZXhwKCkuc3ViKCByb3RhdGVkLnNjYWxlKC0xKS5leHAoKSApLmRpdihpLnNjYWxlKDIpKTtcclxuXHR9XHJcblxyXG5cdGNvcygpIHtcclxuXHRcdC8qXHJcblx0XHRDYWxjdWxhdGUgdGhlIGNvc2luZSBvZiB0aGUgY29tcGxleCBudW1iZXJcclxuXHRcdCovXHJcblx0XHRjb25zdCBpID0gY29tcGxleCgwLCAxKTtcclxuXHRcdGNvbnN0IHJvdGF0ZWQgPSB0aGlzLm11bHQoaSk7XHJcblx0XHRyZXR1cm4gcm90YXRlZC5leHAoKS5hZGQoIHJvdGF0ZWQuc2NhbGUoLTEpLmV4cCgpICkuc2NhbGUoMC41KTtcclxuXHR9XHJcblxyXG5cdHRhbigpIHtcclxuXHRcdC8qXHJcblx0XHRDYWxjdWxhdGUgdGhlIHRhbmdlbnQgb2YgdGhlIGNvbXBsZXggbnVtYmVyXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHRoaXMuc2luKCkuZGl2KHRoaXMuY29zKCkpXHJcblx0fVxyXG5cclxuXHRzaW5oKCkge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgaHlwZXJib2xpYyBzaW5lIG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB0aGlzLmV4cCgpLnN1Yih0aGlzLnNjYWxlKC0xKS5leHAoKSkuc2NhbGUoMC41KTtcclxuXHR9XHJcblxyXG5cdGNvc2goKSB7XHJcblx0XHQvKlxyXG5cdFx0Q2FsY3VsYXRlIHRoZSBoeXBlcmJvbGljIGNvc2luZSBvZiB0aGUgY29tcGxleCBudW1iZXJcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gdGhpcy5leHAoKS5hZGQodGhpcy5zY2FsZSgtMSkuZXhwKCkpLnNjYWxlKDAuNSk7XHJcblx0fVxyXG5cclxuXHR0YW5oKCkge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgaHlwZXJib2xpYyB0YW5nZW50IG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB0aGlzLnNpbmgoKS5kaXYodGhpcy5jb3NoKCkpXHJcblx0fVxyXG5cclxuXHRwb3coeikge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgenRoIHBvd2VyIG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0Ki9cclxuXHRcdGlmICh0aGlzLmVxdWFsc0Vwcyhjb21wbGV4KDAsIDApKSkge1xyXG5cdFx0XHRyZXR1cm4gKHouZXF1YWxzRXBzKGNvbXBsZXgoMSwgMCkpKSA/IGNvbXBsZXgoMSwgMCkgOiBjb21wbGV4KDAsIDApO1xyXG5cdFx0fVxyXG5cclxuXHRcdGNvbnN0IHN1YkFuZyA9IE1hdGguYXRhbjIodGhpcy5pbSwgdGhpcy5yZSk7XHJcblx0XHRjb25zdCBub3JtU3EgPSB0aGlzLm5vcm1TcSgpO1x0XHRcclxuXHRcdGNvbnN0IGFuZyA9IDAuNSAqIHouaW0gKiBNYXRoLmxvZyhub3JtU3EpICsgei5yZSAqIHN1YkFuZztcclxuXHRcdGNvbnN0IG5vcm0gPSBNYXRoLmV4cCgtei5pbSAqIHN1YkFuZykgKiBNYXRoLnBvdyhub3JtU3EsIDAuNSAqIHoucmUpO1xyXG5cdFx0cmV0dXJuIGNvbXBsZXgobm9ybSAqIE1hdGguY29zKGFuZyksIG5vcm0gKiBNYXRoLnNpbihhbmcpKTtcclxuXHR9XHJcblxyXG5cdGNsYW1wKG1pbiwgbWF4KSB7XHJcblx0XHQvKipcclxuXHRcdCAqIGNsYW1wIHRoZSBjb21wbGV4IG51bWJlcidzIG5vcm0gYmV0d2VlbiB0d28gdmFsdWVzXHJcblx0XHQgKi9cclxuXHRcdGNvbnN0IG5vcm0gPSB6Lm5vcm0oKTtcclxuXHRcdGlmICghKG1pbiA8PSBub3JtICYmIG5vcm0gPD0gbWF4KSkge1xyXG5cdFx0XHRyZXR1cm4gdGhpcy5zY2FsZShjbGFtcChub3JtLCBtaW4sIG1heCkgLyBub3JtKTtcclxuXHRcdH1cclxuXHRcdHJldHVybiB0aGlzO1xyXG5cdH1cclxuXHJcblx0ZnJhYygpIHtcclxuXHRcdC8qKlxyXG5cdFx0ICogcmV0dXJuIHRoZSBmcmFjdGlvbmFsIHBhcnQgb2YgZWFjaCBvZiB0aGUgY29tcG9uZW50cyBvZiB0aGUgY29tcGxleCBudW1iZXJcclxuXHRcdCAqL1xyXG5cdFx0cmV0dXJuIGNvbXBsZXgoXHJcblx0XHRcdGZyYWModGhpcy5yZSksXHJcblx0XHRcdGZyYWModGhpcy5pbSksXHJcblx0XHQpO1xyXG5cdH1cclxuXHJcblx0LyogLS0tLS0tLS0tLSBTdGF0aWMgZnVuY3Rpb25zIC0tLS0tLS0tLS0tLS0tLS0tLS0tICovXHJcblxyXG5cdHN0YXRpYyBub3JtKHopIHtcclxuXHRcdHJldHVybiB6Lm5vcm0oKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBub3JtU3Eoeikge1xyXG5cdFx0cmV0dXJuIHoubm9ybVNxKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgaW52KHopIHtcclxuXHRcdHJldHVybiB6LmNvbmooKS5zY2FsZSgxIC8gei5ub3JtU3EoKSk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgc2NhbGUoeiwgaykge1xyXG5cdFx0cmV0dXJuIHouc2NhbGUoayk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgYWRkKHoxLCB6Mikge1xyXG5cdFx0cmV0dXJuIHoxLmFkZCh6Mik7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgc3ViKHoxLCB6Mikge1xyXG5cdFx0cmV0dXJuIHoxLnN1Yih6Mik7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgbXVsdCh6MSwgejIpIHtcclxuXHRcdHJldHVybiB6MS5tdWx0KHoyKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBkaXYoejEsIHoyKSB7XHJcblx0XHRyZXR1cm4gejEuZGl2KHoyKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBwb3coejEsIHoyKSB7XHJcblx0XHRyZXR1cm4gejEucG93KHoyKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBleHAoeikge1xyXG5cdFx0cmV0dXJuIHouZXhwKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgbG4oeikge1xyXG5cdFx0cmV0dXJuIHoubG4oKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBzcXJ0KHopIHtcclxuXHRcdHJldHVybiB6LnNxcnQoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBzaW4oeikge1xyXG5cdFx0LypcclxuXHRcdFJldHVybiB0aGUgc2luZSBvZiB6XHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHouc2luKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgY29zKHopIHtcclxuXHRcdC8qXHJcblx0XHRSZXR1cm4gdGhlIGNvc2luZSBvZiB6XHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHouY29zKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgdGFuKHopIHtcclxuXHRcdC8qXHJcblx0XHRSZXR1cm4gdGhlIHRhbmdlbnQgb2Ygelx0XHRcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gei50YW4oKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBhc2luKHopIHtcclxuXHRcdHJldHVybiB6LmFzaW4oKTtcclxuXHR9IFxyXG5cclxuXHRzdGF0aWMgYWNvcyh6KSB7XHJcblx0XHRyZXR1cm4gei5hY29zKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgYXRhbih6KSB7XHJcblx0XHRyZXR1cm4gei5hdGFuKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgc2luaCh6KSB7XHJcblx0XHQvKlxyXG5cdFx0UmV0dXJuIHRoZSBoeXBlcmJvbGljIHNpbmUgb2YgelxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB6LnNpbmgoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBjb3NoKHopIHtcclxuXHRcdC8qXHJcblx0XHRSZXR1cm4gdGhlIGh5cGVyYm9saWMgY29zaW5lIG9mIHpcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gei5jb3NoKCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgdGFuaCh6KSB7XHJcblx0XHQvKlxyXG5cdFx0UmV0dXJuIHRoZSBoeXBlcmJvbGljIHRhbmdlbnQgb2YgelxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB6LnRhbmgoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBpbmZpbml0ZSh6KSB7XHJcblx0XHQvKipcclxuXHRcdCAqIHJldHVybiB3aGV0aGVyIG9uZSBvciBib3RoIG9mIHRoZSBjb21wb25lbnRzIG9mIHogaXMgaW5maW5pdGUsIG9yIGlmIHRoZSBub3JtIGlzIGluZmluaXRlXHJcblx0XHQgKi9cclxuXHRcdGNvbnN0IG5vcm0gPSB6Lm5vcm0oKTtcclxuXHRcdHJldHVybiAoXHJcblx0XHRcdHoucmUgPT09IEluZmluaXR5IHx8IHoucmUgPT09IC1JbmZpbml0eSBcclxuXHRcdFx0fHwgei5pbSA9PT0gSW5maW5pdHkgfHwgei5pbSA9PT0gLUluZmluaXR5XHJcblx0XHRcdHx8IG5vcm0gPT09IEluZmluaXR5XHJcblx0XHQpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIG5hbih6KSB7XHJcblx0XHQvKipcclxuXHRcdCAqIHJldHVybiB3aGV0aGVyIG9uZSBvciBib3RoIGNvbXBvbmVudHMgb2YgeiBpcyBOYU4gXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIGlzTmFOKHoucmUpIHx8IGlzTmFOKHouaW0pO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGdhbW1hKHopIHtcclxuXHRcdGlmICh6LnJlIDwgMC41KSB7XHJcblx0XHRcdC8vIGdhbW1hKDEteilnYW1tYSh6KSA9IHBpIC8gc2luKHBpKnopXHJcblx0XHRcdHJldHVybiBDb21wbGV4LmRpdiggY29tcGxleChNYXRoLlBJLCAwLjApLCBDb21wbGV4Lm11bHQoei5zY2FsZShNYXRoLlBJKS5zaW4oKSwgQ29tcGxleC5nYW1tYShjb21wbGV4KDEgLSB6LnJlLCAtei5pbSkpKSApO1xyXG5cdFx0fSBlbHNlIHtcclxuXHRcdFx0eiA9IENvbXBsZXguc3ViKHosIGNvbXBsZXgoMSwgMCkpOyAvLyBhY2NvdW50IGZvciBzdHVwaWQgc2hpZnQgYnkgMVxyXG5cdFx0XHRsZXQgeCA9IGNvbXBsZXgocFZhbHVlc1swXSwgMCk7XHJcblx0XHRcdGZvciAobGV0IGk9MTsgaTxwVmFsdWVzLmxlbmd0aDsgaSsrKSB7XHJcblx0XHRcdFx0eCA9IENvbXBsZXguYWRkKHgsIENvbXBsZXguZGl2KCBjb21wbGV4KHBWYWx1ZXNbaV0sIDAuMCksIENvbXBsZXguYWRkKHosIGNvbXBsZXgoaSwgMCkpICkpO1xyXG5cdFx0XHR9XHJcblx0XHRcdGNvbnN0IHQgPSBDb21wbGV4LmFkZCh6LCBjb21wbGV4KDcuNSwgMCkpOyAvLyBnPTcsIGcrMC41XHJcblx0XHRcdHJldHVybiBDb21wbGV4LnBvdyh0LCB6LmFkZChjb21wbGV4KDAuNSwgMCkpKS5tdWx0KHQuc2NhbGUoLTEpLmV4cCgpKS5tdWx0KHgpLnNjYWxlKE1hdGguc3FydCgyICogTWF0aC5QSSkpO1xyXG5cdFx0fVxyXG5cdH1cclxuXHJcblx0c3RhdGljIGJldGEoejEsIHoyKSB7XHJcblx0XHRyZXR1cm4gQ29tcGxleC5kaXYoQ29tcGxleC5tdWx0KENvbXBsZXguZ2FtbWEoejEpLCBDb21wbGV4LmdhbW1hKHoyKSksIENvbXBsZXguZ2FtbWEoejEuYWRkKHoyKSkpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIG1heCh6MSwgejIpIHtcclxuXHRcdHJldHVybiAoejEubm9ybVNxKCkgPCB6Mi5ub3JtU3EoKSkgPyB6MiA6IHoxO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIG1pbih6MSwgejIpIHtcclxuXHRcdHJldHVybiAoejEubm9ybVNxKCkgPiB6Mi5ub3JtU3EoKSkgPyB6MiA6IHoxO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGNsYW1wKHosIG1pbiwgbWF4KSB7XHJcblx0XHRyZXR1cm4gei5jbGFtcChtaW4sIG1heCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgZnJhYygpIHtcclxuXHRcdHJldHVybiB6LmZyYWMoKTtcclxuXHR9XHJcblxyXG5cdC8qIC0tLS0tLS0tLS0tLS0tLSBJbi1wbGFjZSBvcGVyYXRpb25zIC0tLS0tLS0tLS0tLS0tLS0tLS0tLSAqL1xyXG5cclxuXHRpYWRkKHopIHtcclxuXHRcdHRoaXMucmUgKz0gei5yZTtcclxuXHRcdHRoaXMuaW0gKz0gei5pbTtcclxuXHR9XHJcblxyXG5cdGlzdWIoeikge1xyXG5cdFx0dGhpcy5yZSAtPSB6LnJlO1xyXG5cdFx0dGhpcy5pbSAtPSB6LmltO1xyXG5cdH1cclxuXHJcblxyXG59XHJcblxyXG5mdW5jdGlvbiBjb21wbGV4KHJlYWwsIGltYWdpbmFyeSkge1xyXG5cdC8qIGluc3RhbnRpYXRlIGEgQ29tcGxleCB3aXRob3V0IG5ldyBrZXl3b3JkICovXHJcblx0cmV0dXJuIG5ldyBDb21wbGV4KHJlYWwsIGltYWdpbmFyeSk7XHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiBpbnRlZ3JhdGVPdmVyUGFyYW1ldGVyKGYsIGE9MCwgYj0xLCBzYW1wbGVzPTEwMCkge1xyXG5cdC8qXHJcblx0cmV0dXJuIHRoZSBsaW5lIGludGVncmFsIG9mIHRoZSBDb21wbGV4LXZhbHVlZCBmdW5jdGlvbiBmXHJcblx0cGFyYW1ldGVyaXplZCBieSB0aGUgcmVhbCBwYXJhbWV0ZXIgdCBvbiBbYSwgYl0sXHJcblx0dXNpbmcgdGhlIHNwZWNpZmllZCBudW1iZXIgb2Ygc2FtcGxlc1xyXG5cdCovXHJcblx0Y29uc3Qgc3RlcCA9IChiIC0gYSkgLyBzYW1wbGVzO1xyXG5cdGxldCB0ID0gYSArIHN0ZXAgLyAyOyAvLyBtaWRwb2ludCBhcHByb3hpbWF0aW9uXHJcblx0bGV0IHJlc3VsdCA9IGNvbXBsZXgoMCwgMCk7XHJcblx0Zm9yIChsZXQgaT0wOyBpPHNhbXBsZXM7IGkrKykge1xyXG5cdFx0cmVzdWx0LmlhZGQoZih0KS5zY2FsZShzdGVwKSk7XHJcblx0XHR0ICs9IHN0ZXA7XHJcblx0fVxyXG5cdHJldHVybiByZXN1bHQ7XHJcbn1cclxuXHJcbmZ1bmN0aW9uIGludGVncmF0ZU92ZXJDb21wbGV4VmFyaWFibGUoZiwgejAsIHoxLCBzYW1wbGVzPTEwMCkge1xyXG5cdC8qXHJcblx0cmV0dXJuIHRoZSBsaW5lIGludGVncmFsIG9mIHRoZSBDb21wbGV4LXZhbHVlZCBmdW5jdGlvbiBmXHJcblx0b24gdGhlIGxpbmUgYmV0d2VlbiB6MCBhbmQgejEsIHVzaW5nIHRoZSBzcGVjaWZpZWQgbnVtYmVyIG9mIGludGVydmFsc1xyXG5cdCovXHJcblx0Y29uc3Qgc3RlcCA9IENvbXBsZXguc3ViKHoxLCB6MCkuc2NhbGUoMSAvIHNhbXBsZXMpO1xyXG5cdGNvbnN0IHN0ZXBOb3JtID0gc3RlcC5ub3JtKCk7XHJcblx0bGV0IHogPSBDb21wbGV4LmFkZCh6MCwgc3RlcC5zY2FsZSgwLjUpKTtcclxuXHRsZXQgcmVzdWx0ID0gY29tcGxleCgwLCAwKTtcclxuXHRmb3IgKGxldCBpPTA7IGk8c2FtcGxlczsgaSsrKSB7XHJcblx0XHRyZXN1bHQuaWFkZChmKHopLnNjYWxlKHN0ZXBOb3JtKSk7XHJcblx0XHR6LmlhZGQoc3RlcCk7XHJcblx0fVxyXG5cdHJldHVybiByZXN1bHQ7XHJcbn1cclxuXHJcblxyXG5cclxuY2xhc3MgVmVjdG9yIHtcclxuXHJcblx0LypcclxuXHRSZXByZXNlbnRzIGEgdmVjdG9yIGluIENeblxyXG5cdCovXHJcblxyXG5cdGNvbnN0cnVjdG9yKHZhbHVlcykge1xyXG5cdFx0dGhpcy52YWx1ZXMgPSB2YWx1ZXM7XHJcblx0XHR0aGlzLmRpbWVuc2lvbiA9IHZhbHVlcy5sZW5ndGg7XHJcblx0fVxyXG5cclxuXHRnZXQoaW5kZXgpIHtcclxuXHRcdHJldHVybiB0aGlzLnZhbHVlc1tpbmRleF07XHJcblx0fVxyXG5cclxuXHRzY2FsZSh6KSB7XHJcblx0XHQvKlxyXG5cdFx0U2NhbGUgZWFjaCBjb21wb25lbnQgYnkgdGhlIGNvbXBsZXggbnVtYmVyIHpcclxuXHRcdCovXHJcblx0XHRjb25zdCByZXN1bHQgPSBbXTtcclxuXHRcdGZvciAobGV0IGk9MDsgaTx0aGlzLmRpbWVuc2lvbjsgaSsrKSB7XHJcblx0XHRcdHJlc3VsdC5wdXNoKENvbXBsZXgubXVsdCh0aGlzLmdldChpKSwgeikpO1xyXG5cdFx0fVxyXG5cdFx0cmV0dXJuIHZlY3RvcihyZXN1bHQpO1xyXG5cdH1cclxuXHJcblx0cmVhbFNjYWxlKGspIHtcclxuXHRcdC8qXHJcblx0XHRTY2FsZSBlYWNoIGNvbXBvbmVudCBieSB0aGUgcmVhbCBudW1iZXIga1xyXG5cdFx0Ki9cclxuXHRcdGNvbnN0IHJlc3VsdCA9IFtdO1xyXG5cdFx0Zm9yIChsZXQgaT0wOyBpPHRoaXMuZGltZW5zaW9uOyBpKyspIHtcclxuXHRcdFx0cmVzdWx0LnB1c2godGhpcy5nZXQoaSkuc2NhbGUoaykpO1xyXG5cdFx0fVxyXG5cdFx0cmV0dXJuIHZlY3RvcihyZXN1bHQpO1xyXG5cdH1cclxuXHJcblx0YWRkKFopIHtcclxuXHRcdGNvbnN0IHJlc3VsdCA9IFtdO1xyXG5cdFx0Zm9yIChsZXQgaT0wOyBpPHRoaXMuZGltZW5zaW9uOyBpKyspIHtcclxuXHRcdFx0cmVzdWx0LnB1c2goQ29tcGxleC5hZGQodGhpcy5nZXQoaSksIFouZ2V0KGkpKSk7XHJcblx0XHR9XHJcblx0XHRyZXR1cm4gdmVjdG9yKHJlc3VsdCk7XHJcblx0fVxyXG5cdFxyXG5cdGRvdChaKSB7XHJcblx0XHQvKlxyXG5cdFx0Q29tcHV0ZSAoSGVybWl0aWFuKSBkb3QgcHJvZHVjdFxyXG5cdFx0RXF1aXZhbGVudCB0byBzdGFuZGFyZCBkb3QgcHJvZHVjdCBpbiB0aGUgY2FzZSBvZiByZWFsIHZlY3RvcnNcclxuXHRcdCovXHJcblx0XHRjb25zdCByZXN1bHQgPSBjb21wbGV4KDAsIDApO1xyXG5cdFx0Zm9yIChsZXQgaT0wOyBpPHRoaXMuZGltZW5zaW9uOyBpKyspIHtcclxuXHRcdFx0cmVzdWx0LmlhZGQoQ29tcGxleC5tdWx0KHRoaXMuZ2V0KGkpLCBaLmdldChpKS5jb25qKCkpKTtcclxuXHRcdH1cclxuXHRcdHJldHVybiByZXN1bHQ7XHJcblx0fVxyXG5cclxuXHRub3JtKCkge1xyXG5cdFx0LypcclxuXHRcdENvbXB1dGUgbm9ybSBvZiB0aGUgdmVjdG9yXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIE1hdGguc3FydCh0aGlzLmRvdCh0aGlzKS5yZSk7XHJcblx0fVxyXG5cclxuXHRhbmdsZVRvKFopIHtcclxuXHRcdC8qXHJcblx0XHRDb21wdXRlIHRoZSBhbmdsZSBiZXR3ZWVuIHRoZSB2ZWN0b3IgYW5kIHZlY3RvciBaXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIE1hdGguYWNvcyggdGhpcy5kb3QoWikucmUgLyAodGhpcy5ub3JtKCkgKiBaLm5vcm0oKSkgKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBsZXJwKFoxLCBaMiwgdCkge1xyXG5cdFx0LypcclxuXHRcdExpbmVhcmx5IGludGVycG9sYXRlXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIFoxLnJlYWxTY2FsZSgxIC0gdCkuYWRkKFoyLnNjYWxlKHQpKTtcclxuXHR9XHJcblxyXG59XHJcblxyXG5cclxuZnVuY3Rpb24gdmVjdG9yKC4uLnZhbHVlcykge1xyXG5cdC8qXHJcblx0VXRpbGl0eSBmdW5jdGlvbiBmb3IgZWFzeSBpbnN0YW50aWF0aW9uIG9mIFZlY3RvciBvYmplY3RzLlxyXG5cdEFsbG93cyBjb21wb25lbnRzIHRvIGJlIHBhc3NlZCBpbmRpdmlkdWFsbHkgb3IgaW4gYW4gYXJyYXlcclxuXHQqL1xyXG5cdGlmICh2YWx1ZXNbMF0gaW5zdGFuY2VvZiBDb21wbGV4KSB7XHJcblx0XHRyZXR1cm4gbmV3IFZlY3Rvcih2YWx1ZXMpO1xyXG5cdH0gZWxzZSB7XHJcblx0XHRyZXR1cm4gbmV3IFZlY3Rvcih2YWx1ZXNbMF0pO1xyXG5cdH1cclxufVxyXG5cclxuXHJcbm1vZHVsZS5leHBvcnRzID0ge1xyXG5cdENvbXBsZXgsIGNvbXBsZXgsXHJcblx0aW50ZWdyYXRlT3ZlckNvbXBsZXhWYXJpYWJsZSwgaW50ZWdyYXRlT3ZlclBhcmFtZXRlcixcclxuXHRWZWN0b3IsIHZlY3RvcixcclxufTsiXX0=
