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

    "squeeze": "Smoothly compresses a \\( 2\\cdot\\operatorname{length} \\) by \\( 2\\cdot\\operatorname{length} \\) square centered at the origin to the unit circle. The incircle of the square is mapped to a circle of radius \\( 0 \\leq \\operatorname{coverage} \\leq 1 \\).",

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
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi4uLy4uL1VzZXJzL3lhbmtlL0FwcERhdGEvUm9hbWluZy9ucG0vbm9kZV9tb2R1bGVzL2Jyb3dzZXJpZnkvbm9kZV9tb2R1bGVzL2Jyb3dzZXItcGFjay9fcHJlbHVkZS5qcyIsImFwcC9zY29wZS5qcyIsImRvY3MvZ3VpZGVfZ2VuLmpzIiwibWF0aC9jb21wbGV4LmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN0UUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3pJQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24oKXtmdW5jdGlvbiByKGUsbix0KXtmdW5jdGlvbiBvKGksZil7aWYoIW5baV0pe2lmKCFlW2ldKXt2YXIgYz1cImZ1bmN0aW9uXCI9PXR5cGVvZiByZXF1aXJlJiZyZXF1aXJlO2lmKCFmJiZjKXJldHVybiBjKGksITApO2lmKHUpcmV0dXJuIHUoaSwhMCk7dmFyIGE9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitpK1wiJ1wiKTt0aHJvdyBhLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsYX12YXIgcD1uW2ldPXtleHBvcnRzOnt9fTtlW2ldWzBdLmNhbGwocC5leHBvcnRzLGZ1bmN0aW9uKHIpe3ZhciBuPWVbaV1bMV1bcl07cmV0dXJuIG8obnx8cil9LHAscC5leHBvcnRzLHIsZSxuLHQpfXJldHVybiBuW2ldLmV4cG9ydHN9Zm9yKHZhciB1PVwiZnVuY3Rpb25cIj09dHlwZW9mIHJlcXVpcmUmJnJlcXVpcmUsaT0wO2k8dC5sZW5ndGg7aSsrKW8odFtpXSk7cmV0dXJuIG99cmV0dXJuIHJ9KSgpIiwiXHJcbmNvbnN0IHsgY29tcGxleCwgQ29tcGxleCB9ID0gcmVxdWlyZShcIi4uL21hdGgvY29tcGxleC5qc1wiKTtcclxuXHJcblxyXG5jb25zdCBzY29wZSA9IHtcclxuICAgIGJ1aWx0aW46IHtcclxuICAgICAgICBcInpcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwielwiLFxyXG4gICAgICAgICAgICBpc1BhcmFtZXRlcjogdHJ1ZSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIC8vIFwidFwiOiB7XHJcbiAgICAgICAgLy8gICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgIC8vICAgICBzaGFkZXJBbGlhczogXCJ0XCIsXHJcbiAgICAgICAgLy8gICAgIGlzUGFyYW1ldGVyOiB0cnVlLFxyXG4gICAgICAgIC8vIH0sXHJcblxyXG4gICAgICAgIFwibm9ybVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcIm5vcm1DXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibm9ybVNxXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibm9ybVNxQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImFyZ1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImFyZ0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJpbnZcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpbnZDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwiZXhwXCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImV4cENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJsblwiOiB7IFxyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJsbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJzcXJ0XCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInNxcnRDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwic2luXCI6IHsgXHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInNpbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJjb3NcIjogeyBcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY29zQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInRhblwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInRhbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhc2luXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYXNpbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhY29zXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYWNvc0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhdGFuXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiYXRhbkNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJzaW5oXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwic2luaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJjb3NoXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY29zaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJ0YW5oXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwidGFuaENcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJhdGFuaFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImF0YW5oQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInJlXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwicmVDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwiaW1cIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpbUNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJHYW1tYVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcIkdhbW1hQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImJldGFcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJiZXRhQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibWluXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibWluQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibWF4XCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwibWF4Q1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcIndcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwibGVycFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImxlcnBDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwid1wiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSwgXCJ0XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMiB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImNvbmpcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJjb25qQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImNsYW1wXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwiY2xhbXBDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwid1wiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSwgXCJ0XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMiB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImZyYWNcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJmcmFjQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImludmVyc2VTQ1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImludmVyc2VTQ0NcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSwgXCJwXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInNjXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogdHJ1ZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwic2NDXCIsXHJcbiAgICAgICAgICAgIGxvY2FsczogeyBcInpcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAwIH0sIFwicFwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDEgfSB9LFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJwbGFuZVRvUFwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IHRydWUsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInBsYW5lVG9QQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcInBcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgdHlwZTogXCJjb21wbGV4XCIsIGluZGV4OiAxIH0gfSxcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwicFRvUGxhbmVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJwVG9QbGFuZUNcIixcclxuICAgICAgICAgICAgbG9jYWxzOiB7IFwielwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDAgfSwgXCJwXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMSB9IH0sXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInNxdWVlemVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiB0cnVlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJzcXVlZXplQ1wiLFxyXG4gICAgICAgICAgICBsb2NhbHM6IHsgXCJ6XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMCB9LCBcImNvdmVyYWdlXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIHR5cGU6IFwiY29tcGxleFwiLCBpbmRleDogMSB9LCBcImxlbmd0aFwiOiB7IGlzRnVuY3Rpb246IGZhbHNlLCB0eXBlOiBcImNvbXBsZXhcIiwgaW5kZXg6IDIgfSB9LFxyXG4gICAgICAgIH0sXHJcblxyXG4gICAgICAgIFwiaVwiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJpXCIsXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcInBpXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogZmFsc2UsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcInBpQ1wiLFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJ0YXVcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwidHBpQ1wiLFxyXG4gICAgICAgIH0sXHJcbiAgICAgICAgXCJlXCI6IHtcclxuICAgICAgICAgICAgaXNGdW5jdGlvbjogZmFsc2UsXHJcbiAgICAgICAgICAgIHNoYWRlckFsaWFzOiBcImVcIixcclxuICAgICAgICB9LFxyXG4gICAgICAgIFwicmVhbEJvdW5kc1wiOiB7XHJcbiAgICAgICAgICAgIGlzRnVuY3Rpb246IGZhbHNlLFxyXG4gICAgICAgICAgICBzaGFkZXJBbGlhczogXCJ4Qm91bmRzXCIsXHJcbiAgICAgICAgfSxcclxuICAgICAgICBcImltYWdCb3VuZHNcIjoge1xyXG4gICAgICAgICAgICBpc0Z1bmN0aW9uOiBmYWxzZSxcclxuICAgICAgICAgICAgc2hhZGVyQWxpYXM6IFwieUJvdW5kc1wiLFxyXG4gICAgICAgIH0sXHJcblxyXG4gICAgICAgIC8qIGRhdGF0eXBlcyAoZm9yIHR5cGUgYW5ub3RhdGlvbikgKi9cclxuICAgICAgICBcImNvbXBsZXhcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgaXNUeXBlOiB0cnVlLCB9LFxyXG4gICAgICAgIFwiYXJyYXlcIjogeyBpc0Z1bmN0aW9uOiBmYWxzZSwgaXNUeXBlOiB0cnVlLCB9LFxyXG4gICAgICAgIFwibWF0cml4XCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIGlzVHlwZTogdHJ1ZSwgfSxcclxuICAgICAgICBcImZ1bmN0aW9uXCI6IHsgaXNGdW5jdGlvbjogZmFsc2UsIGlzVHlwZTogdHJ1ZSwgfSxcclxuICAgIH0sXHJcblxyXG4gICAgdXNlckdsb2JhbDoge1xyXG4gICAgfSxcclxufVxyXG5cclxuXHJcbmNvbnN0IGRlZmF1bHRWYWx1ZVNjb3BlID0ge1xyXG4gICAgXCJ6XCI6IGNvbXBsZXgoMSwgMCksXHJcbiAgICBcIm5vcm1cIjogKHopID0+IGNvbXBsZXgoei5ub3JtKCksIDApLFxyXG4gICAgXCJub3JtU3FcIjogKHopID0+IGNvbXBsZXgoei5ub3JtU3EoKSwgMCksXHJcbiAgICBcImFyZ1wiOiAoeikgPT4gY29tcGxleCh6LmFyZygpLCAwKSxcclxuICAgIFwiaW52XCI6IENvbXBsZXguaW52LFxyXG4gICAgXCJleHBcIjogQ29tcGxleC5leHAsXHJcbiAgICBcImxuXCI6IENvbXBsZXgubG4sXHJcbiAgICBcInNxcnRcIjogQ29tcGxleC5zcXJ0LFxyXG4gICAgXCJzaW5cIjogQ29tcGxleC5zaW4sXHJcbiAgICBcImNvc1wiOiBDb21wbGV4LmNvcyxcclxuICAgIFwidGFuXCI6IENvbXBsZXgudGFuLFxyXG4gICAgXCJhY29zXCI6IENvbXBsZXguYWNvcyxcclxuICAgIFwiYXNpblwiOiBDb21wbGV4LmFzaW4sXHJcbiAgICBcImF0YW5cIjogQ29tcGxleC5hdGFuLFxyXG4gICAgXCJzaW5oXCI6IENvbXBsZXguc2luaCxcclxuICAgIFwiY29zaFwiOiBDb21wbGV4LmNvc2gsXHJcbiAgICBcInRhbmhcIjogQ29tcGxleC50YW5oLFxyXG4gICAgXCJhdGFuaFwiOiAoeikgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInJlXCI6ICh6KSA9PiBjb21wbGV4KHoucmUsIDApLFxyXG4gICAgXCJpbVwiOiAoeikgPT4gY29tcGxleCh6LmltLCAwKSxcclxuICAgIFwiR2FtbWFcIjogQ29tcGxleC5nYW1tYSxcclxuICAgIFwiYmV0YVwiOiBDb21wbGV4LmJldGEsXHJcbiAgICBcIm1pblwiOiBDb21wbGV4Lm1pbixcclxuICAgIFwibWF4XCI6IENvbXBsZXgubWF4LFxyXG4gICAgXCJsZXJwXCI6ICh6MSwgejIsIHQpID0+IENvbXBsZXgubXVsdChjb21wbGV4KDEsIDApLnN1Yih0KSwgejEpLmFkZCh6Mi5tdWx0KHQpKSxcclxuICAgIFwiY29ualwiOiAoeikgPT4gei5jb25qKCksXHJcbiAgICBcImNsYW1wXCI6ICh6LCBtaW4sIG1heCkgPT4gQ29tcGxleC5jbGFtcCh6LCBtaW4ubm9ybSgpLCBtYXgubm9ybSgpKSxcclxuICAgIFwiZnJhY1wiOiBDb21wbGV4LmZyYWMsXHJcbiAgICBcImludmVyc2VTQ1wiOiAoeiwgcCkgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInNjXCI6ICh6LCBwKSA9PiBjb21wbGV4KDEsIDApLCAvLyBub3QgaW1wbGVtZW50ZWQgZm9yIHA1IG1vZGVcclxuICAgIFwicGxhbmVUb1BcIjogKHosIHApID0+IGNvbXBsZXgoMSwgMCksIC8vIG5vdCBpbXBsZW1lbnRlZCBmb3IgcDUgbW9kZVxyXG4gICAgXCJwVG9QbGFuZVwiOiAoeiwgcCkgPT4gY29tcGxleCgxLCAwKSwgLy8gbm90IGltcGxlbWVudGVkIGZvciBwNSBtb2RlXHJcbiAgICBcInNxdWVlemVcIjogKHosIGNvdmVyYWdlLCBsZW5ndGgpID0+IGNvbXBsZXgoMSwgMCksIC8vIG5vdCBpbXBsZW1lbnRlZCBmb3IgcDUgbW9kZVxyXG4gIFxyXG4gICAgXCJpXCI6IGNvbXBsZXgoMCwgMSksXHJcbiAgICBcInBpXCI6IGNvbXBsZXgoTWF0aC5QSSwgMCksXHJcbiAgICBcInRhdVwiOiBjb21wbGV4KDIgKiBNYXRoLlBJLCAwKSxcclxuICAgIFwiZVwiOiBjb21wbGV4KE1hdGguRSwgMCksXHJcbiAgICBcInJlYWxCb3VuZHNcIjogY29tcGxleCgwLCAwKSxcclxuICAgIFwiaW1hZ0JvdW5kc1wiOiBjb21wbGV4KDAsIDApLFxyXG59O1xyXG5cclxubGV0IHZhbHVlU2NvcGUgPSB7fTtcclxuXHJcbm1vZHVsZS5leHBvcnRzID0ge1xyXG4gICAgc2NvcGUsIGRlZmF1bHRWYWx1ZVNjb3BlLCB2YWx1ZVNjb3BlLFxyXG59OyIsIlxyXG5cclxuXHJcbmNvbnN0IHNjb3BlcyA9IHJlcXVpcmUoXCIuLi9hcHAvc2NvcGUuanNcIik7XHJcbmNvbnN0IHNjb3BlID0gc2NvcGVzLnNjb3BlLmJ1aWx0aW47XHJcblxyXG5jb25zdCBkZXNjcmlwdGlvbnMgPSB7XHJcbiAgICBcIkdhbW1hXCI6IFwiVGhlIEdhbW1hIGZ1bmN0aW9uIFxcXFwoIFxcXFxHYW1tYSh6KSBcXFxcKS5cIixcclxuXHJcbiAgICBcImFjb3NcIjogXCJUaGUgaW52ZXJzZSBjb3NpbmUgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJhcmdcIjogXCJUaGUgYXJndW1lbnQgb3IgYW5nbGUgb2YgYSBjb21wbGV4IG51bWJlci5cIixcclxuXHJcbiAgICBcImFycmF5XCI6IFwiVGhlIGFycmF5IGRhdGF0eXBlLiBUaGlzIGlzIGN1cnJlbnRseSBub3QgeWV0IGltcGxlbWVudGVkLCBidXQgdGhlIGtleXdvcmQgaXMgYWxyZWFkeSByZXNlcnZlZCBpbiB0aGUgY2FsY3VsYXRvci5cIixcclxuXHJcbiAgICBcImFzaW5cIjogXCJUaGUgaW52ZXJzZSBzaW5lIGZ1bmN0aW9uLlwiLFxyXG5cclxuICAgIFwiYXRhblwiOiBcIlRoZSBpbnZlcnNlIHRhbmdlbnQgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJhdGFuaFwiOiBcIlRoZSBpbnZlcnNlIGh5cGVyYm9saWMgdGFuZ2VudCBmdW5jdGlvbi5cIixcclxuXHJcbiAgICBcImJldGFcIjogXCJUaGUgYmV0YSBmdW5jdGlvbiBcXFxcKCBcXFxcYmV0YSh6LCB3KSBcXFxcKS5cIixcclxuXHJcbiAgICBcImNsYW1wXCI6IFwiTWFrZSBzdXJlIHRoYXQgXFxcXCggfHd8IFxcXFxsZXEgfHp8IFxcXFxsZXEgfHR8ICBcXFxcKS4gSWYgc28sIHJldHVybiBcXFxcKCB6IFxcXFwpLiBJZiBub3QsIHJldHVybiBcXFxcKCB6IFxcXFwpIHNjYWxlZCB0byB0aGUgbm9ybSBvZiB0aGUgY2xvc2VyIGJvdW5kLiBGb3IgZXhhbXBsZSwgXFxcXCggXFxcXG9wZXJhdG9ybmFtZXtjbGFtcH0oMS4xLCAwLCAxKT0xIFxcXFwpLlwiLFxyXG5cclxuICAgIFwiY29tcGxleFwiOiBcIlRoZSBjb21wbGV4IGRhdGF0eXBlLiBJdCByZXByZXNlbnRzIGEgc2luZ2xlIG51bWJlciBvZiB0aGUgZm9ybSBcXFxcKCBhK2liIFxcXFwpLlwiLFxyXG5cclxuICAgIFwiY29ualwiOiBcIlRoZSBjb21wbGV4IGNvbmp1Z2F0ZSBvZiBhIGNvbXBsZXggbnVtYmVyLCB0aGF0IGlzLCB0aGUgcmVmbGVjdGlvbiBvZiB0aGUgbnVtYmVyIGFjcm9zcyB0aGUgcmVhbCBheGlzLlwiLFxyXG5cclxuICAgIFwiY29zXCI6IFwiVGhlIGNvc2luZSBmdW5jdGlvbi5cIixcclxuXHJcbiAgICBcImNvc2hcIjogXCJUaGUgaHlwZXJib2xpYyBjb3NpbmUgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJlXCI6IFwiRXVsZXIncyBjb25zdGFudCwgdGhlIGJhc2Ugb2YgdGhlIG5hdHVyYWwgbG9nYXJpdGhtLCBcXFxcKCBlIFxcXFxhcHByb3ggMi43MTggXFxcXCkuXCIsXHJcblxyXG4gICAgXCJleHBcIjogXCJUaGUgZXhwb25lbnRpYWwgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJmcmFjXCI6IFwiUmV0dXJucyB0aGUgZnJhY3Rpb25hbCBwYXJ0ICh0aGUgcGFydCBhZnRlciB0aGUgZGVjaW1hbCBwb2ludCkgb2YgYm90aCB0aGUgcmVhbCBhbmQgaW1hZ2luYXJ5IGNvbXBvbmVudHMgb2YgYSBjb21wbGV4IG51bWJlci4gRm9yIGV4YW1wbGUsIFxcXFwoIFxcXFxvcGVyYXRvcm5hbWV7ZnJhY30oMS4xMjMgKyA0LjU2N2kpPTAuMTIzKzAuNTY3aSBcXFxcKS5cIixcclxuXHJcbiAgICBcImZ1bmN0aW9uXCI6IFwiVGhlIGZ1bmN0aW9uIGRhdGF0eXBlXCIsXHJcbiAgICBcclxuICAgIFwiaVwiOiBcIlRoZSBpbWFnaW5hcnkgdW5pdCBcXFxcKCBpIDo9IHggXFxcXGluIFxcXFxtYXRoYmJ7Un1cXFxcbGVmdFt4XFxcXHJpZ2h0XS8oeF4yKzEpIFxcXFwpLlwiLFxyXG5cclxuICAgIFwiaW1cIjogXCJUaGUgaW1hZ2luYXJ5IHBhcnQgb2YgYSBjb21wbGV4IG51bWJlci5cIixcclxuXHJcbiAgICBcImltYWdCb3VuZHNcIjogXCJUaGUgY3VycmVudCB5LWF4aXMgYm91bmRzIGluIHBsYW5lIG1vZGUuIFRoZSBmb3JtYXQgaXMgXFxcXCggXFxcXG9wZXJhdG9ybmFtZXtpbWFnQm91bmRzfSA9IFxcXFxvcGVyYXRvcm5hbWV7eU1pbn0gKyBcXFxcb3BlcmF0b3JuYW1le3lNYXh9XFxcXGNkb3QgaSBcXFxcKS5cIixcclxuXHJcbiAgICBcImludlwiOiBcIlRoZSBtdWx0aXBsaWNhdGl2ZSBpbnZlcnNlIG9mIGEgY29tcGxleCBudW1iZXIuXCIsXHJcblxyXG4gICAgXCJpbnZlcnNlU0NcIjogXCJDb25zdHJ1Y3Rpb24gb2YgdGhlIGludmVyc2UgU2Nod2Fyei1DaHJpc3RvZmZlbCBtYXAgZm9yIGEgcmVndWxhciBcXFxcKCBwIFxcXFwpLWdvbi4gVGhpcyBjb25mb3JtYWxseSBtYXBzIHRoZSB1bml0IGRpc2sgdG8gdGhlIFxcXFwoIHAgXFxcXCktZ29uLlwiLFxyXG5cclxuICAgIFwibGVycFwiOiBcIkxpbmVhcmx5IGludGVycG9sYXRlcyBcXFxcKCB6IFxcXFwpIGJldHdlZW4gXFxcXCggdyBcXFxcKSBhbmQgXFxcXCggdCBcXFxcKS5cIixcclxuXHJcbiAgICBcImxuXCI6IFwiVGhlIG5hdHVyYWwgKGJhc2UgXFxcXCggZSBcXFxcKSBsb2dhcml0aG0uXCIsXHJcblxyXG4gICAgXCJtYXRyaXhcIjogXCJUaGUgbWF0cml4IGRhdGF0eXBlLiBUaGlzIGlzIGN1cnJlbnRseSBub3QgeWV0IGltcGxlbWVudGVkLCBidXQgdGhlIGtleXdvcmQgaXMgYWxyZWFkeSByZXNlcnZlZCBpbiB0aGUgY2FsY3VsYXRvci5cIixcclxuICAgIFxyXG4gICAgXCJtYXhcIjogXCJSZXR1cm5zIHRoZSBhcmd1bWVudCB3aXRoIHRoZSBncmVhdGVyIG5vcm0uXCIsXHJcblxyXG4gICAgXCJtaW5cIjogXCJSZXR1cm5zIHRoZSBhcmd1bWVudCB3aXRoIHRoZSBsZXNzZXIgbm9ybS5cIixcclxuICAgIFxyXG4gICAgXCJub3JtXCI6IFwiVGhlIEV1Y2xpZGVhbiBub3JtIChtYWduaXR1ZGUpIG9mIGEgY29tcGxleCBudW1iZXIuXCIsXHJcblxyXG4gICAgXCJub3JtU3FcIjogXCJUaGUgc3F1YXJlIG9mIHRoZSBFdWNsaWRlYW4gbm9ybSAobWFnbml0dWRlKSBvZiBhIGNvbXBsZXggbnVtYmVyLlwiLFxyXG5cclxuICAgIFwicFRvUGxhbmVcIjogXCJBIGNvbmZvcm1hbCBtYXAgZnJvbSBhIHJlZ3VsYXIgdW5pdC1yYWRpdXMgXFxcXCggcCBcXFxcKS1nb24gdG8gdGhlIGVudGlyZXR5IG9mIHRoZSBjb21wbGV4IHBsYW5lLlwiLFxyXG5cclxuICAgIFwicGlcIjogXCJUaGUgbWF0aGVtYXRpY2FsIGhhbGYtY2lyY2xlIGNvbnN0YW50IFxcXFwoIFxcXFxwaSBcXFxcYXBwcm94IDMuMTQyIFxcXFwpLlwiLFxyXG5cclxuICAgIFwicGxhbmVUb1BcIjogXCJBIGNvbmZvcm1hbCBtYXAgZnJvbSB0aGUgY29tcGxleCBwbGFuZSB0byBhIHJlZ3VsYXIgdW5pdC1yYWRpdXMgXFxcXCggcCBcXFxcKS1nb24uXCIsXHJcblxyXG4gICAgXCJyZVwiOiBcIlRoZSByZWFsIHBhcnQgb2YgYSBjb21wbGV4IG51bWJlci5cIixcclxuXHJcbiAgICBcInJlYWxCb3VuZHNcIjogXCJUaGUgY3VycmVudCB4LWF4aXMgYm91bmRzIGluIHBsYW5lIG1vZGUuIFRoZSBmb3JtYXQgaXMgXFxcXCggXFxcXG9wZXJhdG9ybmFtZXtpbWFnQm91bmRzfSA9IFxcXFxvcGVyYXRvcm5hbWV7eE1pbn0gKyBcXFxcb3BlcmF0b3JuYW1le3hNYXh9XFxcXGNkb3QgaSBcXFxcKS5cIixcclxuXHJcbiAgICBcInNjXCI6IFwiQ29uc3RydWN0aW9uIG9mIHRoZSBTY2h3YXJ6LUNocmlzdG9mZmVsIG1hcCBmb3IgYSByZWd1bGFyIFxcXFwoIHAgXFxcXCktZ29uLiBUaGlzIGNvbmZvcm1sbHkgbWFwcyBhIHJlZ3VsYXIgXFxcXCggcCBcXFxcKS1nb24gdG8gdGhlIHVuaXQgZGlzay5cIixcclxuXHJcbiAgICBcInNpblwiOiBcIlRoZSBzaW5lIGZ1bmN0aW9uLlwiLFxyXG5cclxuICAgIFwic2luaFwiOiBcIlRoZSBoeXBlcmJvbGljIHNpbmUgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJzcXJ0XCI6IFwiVGhlIHNxdWFyZSByb290IGZ1bmN0aW9uLlwiLFxyXG5cclxuICAgIFwic3F1ZWV6ZVwiOiBcIlNtb290aGx5IGNvbXByZXNzZXMgYSBcXFxcKCAyXFxcXGNkb3RcXFxcb3BlcmF0b3JuYW1le2xlbmd0aH0gXFxcXCkgYnkgXFxcXCggMlxcXFxjZG90XFxcXG9wZXJhdG9ybmFtZXtsZW5ndGh9IFxcXFwpIHNxdWFyZSBjZW50ZXJlZCBhdCB0aGUgb3JpZ2luIHRvIHRoZSB1bml0IGNpcmNsZS4gVGhlIGluY2lyY2xlIG9mIHRoZSBzcXVhcmUgaXMgbWFwcGVkIHRvIGEgY2lyY2xlIG9mIHJhZGl1cyBcXFxcKCAwIFxcXFxsZXEgXFxcXG9wZXJhdG9ybmFtZXtjb3ZlcmFnZX0gXFxcXGxlcSAxIFxcXFwpLlwiLFxyXG5cclxuICAgIFwidGFuXCI6IFwiVGhlIHRhbmdlbnQgZnVuY3Rpb24uXCIsXHJcblxyXG4gICAgXCJ0YW5oXCI6IFwiVGhlIGh5cGVyYm9saWMgdGFuZ2VudCBmdW5jdGlvbi5cIixcclxuXHJcbiAgICBcInRhdVwiOiBcIlRoZSBtYXRoZW1hdGljYWwgY2lyY2xlIGNvbnN0YW50IFxcXFwoIFxcXFx0YXUgXFxcXGFwcHJveCA2LjI4MyBcXFxcKS5cIixcclxuXHJcbiAgICBcInpcIjogXCJBIHBhcmFtZXRlciByZXByZXNlbnRpbmcgYSBjb21wbGV4IG51bWJlciB0aGF0IHJhbmdlcyBhY3Jvc3MgdGhlIGVudGlyZSBwbGFuZS5cIixcclxufTtcclxuXHJcblxyXG5jb25zdCBnZW5lcmF0ZURlc2NyaXB0aW9ucyA9ICgpID0+IHtcclxuICAgIGNvbnN0IGtleXMgPSBPYmplY3Qua2V5cyhzY29wZSkuc29ydCgpO1xyXG4gICAgbGV0IHJlc3VsdCA9IFwiXCI7XHJcblxyXG4gICAgZm9yIChjb25zdCBrZXkgb2Yga2V5cykge1xyXG4gICAgICAgIGNvbnN0IGRlc2NyaXB0aW9uID0gZGVzY3JpcHRpb25zW2tleV0gPz8gXCJObyBkZXNjcmlwdGlvbiBnaXZlblwiO1xyXG4gICAgICAgIGxldCBhcmd1bWVudHMgPSBcIlwiO1xyXG4gICAgICAgIGlmIChzY29wZVtrZXldLmlzRnVuY3Rpb24pIHtcclxuICAgICAgICAgICAgYXJndW1lbnRzID0gXCJBcmd1bWVudHM6IFwiO1xyXG4gICAgICAgICAgICBjb25zdCBsb2NhbHMgPSBPYmplY3Qua2V5cyhzY29wZVtrZXldLmxvY2Fscykuc29ydChsb2NhbCA9PiBsb2NhbC5pbmRleCk7XHJcbiAgICAgICAgICAgIGlmIChsb2NhbHMubGVuZ3RoID09PSAwKSBhcmd1bWVudHMgKz0gXCJub25lXCI7XHJcbiAgICAgICAgICAgIGZvciAobGV0IGk9MDsgaTxsb2NhbHMubGVuZ3RoOyBpKyspIHtcclxuICAgICAgICAgICAgICAgIGNvbnN0IGxvY2FsID0gbG9jYWxzW2ldO1xyXG4gICAgICAgICAgICAgICAgYXJndW1lbnRzICs9IGxvY2FsICsgYCAoJHtzY29wZVtrZXldLmxvY2Fsc1tsb2NhbF0udHlwZX0pYDtcclxuICAgICAgICAgICAgICAgIGlmIChpICE9PSBsb2NhbHMubGVuZ3RoIC0gMSkgYXJndW1lbnRzICs9IFwiLCBcIjtcclxuICAgICAgICAgICAgfVxyXG4gICAgICAgICAgICBhcmd1bWVudHMgKz0gXCI8YnI+XCI7XHJcbiAgICAgICAgfVxyXG4gICAgICAgIHJlc3VsdCArPSBgPGRpdiBjbGFzcz1cImRlc2NyaXB0aW9uLWVudHJ5XCI+PHNwYW4gc3R5bGU9XCJmb250LXdlaWdodDogYm9sZFwiPiR7a2V5fTwvc3Bhbj48YnI+JHthcmd1bWVudHN9RGVzY3JpcHRpb246ICR7ZGVzY3JpcHRpb259PC9kaXY+YDtcclxuICAgIH1cclxuXHJcbiAgICByZXR1cm4gcmVzdWx0O1xyXG59O1xyXG5cclxuY29uc3Qgb3V0cHV0RGVzY3JpcHRpb25zID0gKHRhcmdldElkKSA9PiB7XHJcbiAgICBjb25zdCB0YXJnZXRFbGVtZW50ID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihgIyR7dGFyZ2V0SWR9YCk7XHJcbiAgICB0YXJnZXRFbGVtZW50LmlubmVySFRNTCA9IGdlbmVyYXRlRGVzY3JpcHRpb25zKCk7XHJcbn1cclxuXHJcbmNvbnN0IHRvZ2dsZURlc2NyaXB0aW9ucyA9ICgpID0+IHtcclxuICAgIGNvbnN0IGRlc2NEaXYgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yKFwiI2Rlc2NyaXB0aW9uLWNvbnRhaW5lclwiKTtcclxuICAgIGNvbnN0IGNvbGxhcHNlQnRuID0gZG9jdW1lbnQucXVlcnlTZWxlY3RvcihcIiNjb2xsYXBzZS1idG5cIik7XHJcbiAgICBpZiAoZGVzY0Rpdi5zdHlsZS5kaXNwbGF5ID09PSBcIm5vbmVcIikge1xyXG4gICAgICAgIGRlc2NEaXYuc3R5bGUuZGlzcGxheSA9IFwiYmxvY2tcIjtcclxuICAgICAgICBjb2xsYXBzZUJ0bi5pbm5lclRleHQgPSBcIkNvbGxhcHNlXCI7XHJcbiAgICB9IGVsc2Uge1xyXG4gICAgICAgIGRlc2NEaXYuc3R5bGUuZGlzcGxheSA9IFwibm9uZVwiO1xyXG4gICAgICAgIGNvbGxhcHNlQnRuLmlubmVyVGV4dCA9IFwiRXhwYW5kXCI7XHJcbiAgICB9XHJcbn1cclxuXHJcbm91dHB1dERlc2NyaXB0aW9ucyhcImRlc2NyaXB0aW9uLWNvbnRhaW5lclwiKTtcclxud2luZG93LnRvZ2dsZURlc2NyaXB0aW9ucyA9IHRvZ2dsZURlc2NyaXB0aW9uczsiLCJjb25zdCBFUFNJTE9OID0gMC4wMDAwMDE7XHJcblxyXG5jb25zdCBwVmFsdWVzID0gW1xyXG5cdDAuOTk5OTk5OTk5OTk5ODA5OTMsXHJcblx0Njc2LjUyMDM2ODEyMTg4NTEsXHJcblx0LTEyNTkuMTM5MjE2NzIyNDAyOCxcclxuXHQ3NzEuMzIzNDI4Nzc3NjUzMTMsXHJcblx0LTE3Ni42MTUwMjkxNjIxNDA1OSxcclxuXHQxMi41MDczNDMyNzg2ODY5MDUsXHJcblx0LTAuMTM4NTcxMDk1MjY1NzIwMTIsXHJcblx0OS45ODQzNjk1NzgwMTk1NzE2ZS02LFxyXG5cdDEuNTA1NjMyNzM1MTQ5MzExNmUtN1xyXG5dO1xyXG5cclxuXHJcbmNvbnN0IGZyYWMgPSAoeCkgPT4ge1xyXG5cdHJldHVybiB4IC0gKHggPiAwID8gTWF0aC5mbG9vcih4KSA6IE1hdGguY2VpbCh4KSk7XHJcbn07XHJcblxyXG5jb25zdCBjbGFtcCA9ICh4LCBtaW4sIG1heCkgPT4ge1xyXG5cdHJldHVybiBNYXRoLm1pbihtYXgsIE1hdGgubWF4KG1pbiwgeCkpO1xyXG59O1xyXG5cclxuXHJcbmNsYXNzIENvbXBsZXgge1xyXG5cclxuXHQvKlxyXG5cdENsYXNzIGZvciByZXByZXNlbnRpbmcgY29tcGxleCBudW1iZXJzIG9mIHRoZSBmb3JtIGEgKyBiaVxyXG5cdCovXHJcblxyXG5cdGNvbnN0cnVjdG9yKHJlYWwsIGltYWdpbmFyeSkge1xyXG5cdFx0dGhpcy5yZSA9IHJlYWw7XHJcblx0XHR0aGlzLmltID0gaW1hZ2luYXJ5O1xyXG5cdH1cclxuXHJcblx0Y29uaigpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBjb21wbGV4IGNvbmp1Z2F0ZSAqL1xyXG5cdFx0cmV0dXJuIG5ldyBDb21wbGV4KHRoaXMucmUsIC10aGlzLmltKTtcclxuXHR9XHJcblxyXG5cdG5vcm0oKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgbm9ybSAobW9kdWx1cyksIGFzIGEgcmVhbCBudW1iZXIgKi9cclxuXHRcdHJldHVybiBNYXRoLnNxcnQodGhpcy5yZSAqIHRoaXMucmUgKyB0aGlzLmltICogdGhpcy5pbSk7XHJcblx0fVxyXG5cclxuXHRub3JtU3EoKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgc3F1YXJlIG9mIHRoZSBub3JtIChtb2R1bHVzKSwgYXMgYSByZWFsIG51bWJlciAqL1xyXG5cdFx0cmV0dXJuIHRoaXMucmUgKiB0aGlzLnJlICsgdGhpcy5pbSAqIHRoaXMuaW07XHJcblx0fVxyXG5cclxuXHRhcmcoKSB7XHJcblx0XHQvKlxyXG5cdFx0Q29tcHV0ZXMgdGhlIGFuZ2xlIChhcmd1bWVudCksIGFzIGEgcmVhbCBudW1iZXIgbWVhc3VyZWQgaW4gcmFkaWFuc1xyXG5cdFx0MCA8PSBhcmcoeikgPCAyICogcGlcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gKE1hdGguYXRhbjIodGhpcy5pbSwgdGhpcy5yZSkgKyAyICogTWF0aC5QSSkgJSAoMiAqIE1hdGguUEkpO1xyXG5cdH1cclxuXHJcblx0dW5pdCgpIHtcclxuXHRcdC8qIENvbXB1dGVzIGEgdW5pdCBtb2R1bHVzIGNvbXBsZXggbnVtYmVyIGluIHRoZSBkaXJlY3Rpb24gb2YgdGhpcyBjb21wbGV4IG51bWJlciAqL1xyXG5cdFx0cmV0dXJuIHRoaXMuc2NhbGUoMSAvIHRoaXMubm9ybSgpKTtcclxuXHR9XHJcblxyXG5cdHNjYWxlKGspIHtcclxuXHRcdC8qIFNjYWxlcyBlYWNoIGNvbXBvbmVudCBieSB0aGUgcmVhbCBjb25zdGFudCBrICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSAqIGssIHRoaXMuaW0gKiBrKTtcclxuXHR9XHJcblxyXG5cdGFkZCh6KSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgc3VtIG9mIHRoaXMgY29tcGxleCBudW1iZXIgYW5kIHogKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleCh0aGlzLnJlICsgei5yZSwgdGhpcy5pbSArIHouaW0pO1xyXG5cdH1cclxuXHJcblx0c3ViKHopIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBkaWZmZXJlbmNlIG9mIHRoaXMgY29tcGxleCBudW1iZXIgYW5kIHogKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleCh0aGlzLnJlIC0gei5yZSwgdGhpcy5pbSAtIHouaW0pO1x0XHJcblx0fVxyXG5cclxuXHRtdWx0KHopIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBwcm9kdWN0IG9mIHRoaXMgY29tcGxleCBudW1iZXIgYW5kIHogKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleCh0aGlzLnJlICogei5yZSAtIHRoaXMuaW0gKiB6LmltLCB0aGlzLnJlICogei5pbSArIHRoaXMuaW0gKiB6LnJlKTtcclxuXHR9XHJcblxyXG5cdGVNdWx0KHopIHtcclxuXHRcdC8qIGVsZW1lbnR3aXNlIG11bHRpcGxpY2F0aW9uIG9mIHRoaXMgY29tcGxleCBudW1iZXIgYW5kIHogKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleCh0aGlzLnJlICogei5yZSwgdGhpcy5pbSAqIHouaW0pO1xyXG5cdH1cclxuXHJcblx0aW52KCkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIHJlY2lwcm9jYWwgKGludmVyc2UpICovXHJcblx0XHRyZXR1cm4gdGhpcy5jb25qKCkuc2NhbGUoMSAvIHRoaXMubm9ybVNxKCkpO1xyXG5cdH1cclxuXHJcblx0ZGl2KHopIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBxdW90aWVudCBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gdGhpcy5tdWx0KHouaW52KCkpO1xyXG5cdH1cclxuXHJcblx0cGVycCgpIHtcclxuXHRcdC8qIENvbXB1dGVzIGFuIG9ydGhvZ29uYWwgY29tcGxleCBudW1iZXIgb2YgdGhlIHNhbWUgbWFnbml0dWRlICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgoLXRoaXMuaW0sIHRoaXMucmUpO1xyXG5cdH1cclxuXHJcblx0c3FydCgpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBwcmluY2lwYWwgYnJhbmNoIG9mIHRoZSBzcXVhcmUgcm9vdCAqL1xyXG5cdFx0Y29uc3Qgbm9ybVNxcnQgPSBNYXRoLnNxcnQodGhpcy5ub3JtKCkpO1xyXG5cdFx0Y29uc3QgaGFsZkFyZyA9IDAuNSAqIHRoaXMuYXJnKCk7XHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgobm9ybVNxcnQgKiBNYXRoLmNvcyhoYWxmQXJnKSwgbm9ybVNxcnQgKiBNYXRoLnNpbihoYWxmQXJnKSk7XHJcblx0fVxyXG5cclxuXHRzcXVhcmUoKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgc3F1YXJlICovXHJcblx0XHRyZXR1cm4gbmV3IENvbXBsZXgodGhpcy5yZSAqIHRoaXMucmUgLSB0aGlzLmltICogdGhpcy5pbSwgMiAqIHRoaXMucmUgKiB0aGlzLmltKTtcclxuXHR9XHJcblxyXG5cdGV4cCgpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBleHBvbmVudGlhbCBmdW5jdGlvbiBvZiB0aGlzIGNvbXBsZXggbnVtYmVyICovXHJcblx0XHRjb25zdCBtYWcgPSBNYXRoLmV4cCh0aGlzLnJlKTtcclxuXHRcdHJldHVybiBuZXcgQ29tcGxleChtYWcgKiBNYXRoLmNvcyh0aGlzLmltKSwgbWFnICogTWF0aC5zaW4odGhpcy5pbSkpO1xyXG5cdH1cclxuXHJcblx0bG4oKSB7XHJcblx0XHQvKiBDb21wdXRlcyB0aGUgcHJpbmNpcGFsIGJyYW5jaCBvZiB0aGUgbmF0dXJhbCBsb2cgKi9cclxuXHRcdHJldHVybiBuZXcgQ29tcGxleChNYXRoLmxvZyh0aGlzLm5vcm0oKSksIHRoaXMuYXJnKCkpO1xyXG5cdH1cclxuXHJcblx0YWNvcygpIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBwcmluY2lwYWwgYnJhbmNoIG9mIHRoZSBpbnZlcnNlIGNvc2luZSAqL1xyXG5cdFx0cmV0dXJuIHRoaXMuYWRkKHRoaXMuc3F1YXJlKCkuc3ViKG5ldyBDb21wbGV4KDEsIDApKS5zcXJ0KCkpLmxuKCkuZGl2KGNvbXBsZXgoMCwgMSkpO1xyXG5cdH1cclxuXHJcblx0YXNpbigpIHtcclxuXHRcdC8qKiBjb21wdXRlcyB0aGUgcHJpbmNpcGxlIGJyYW5jaCBvZiB0aGUgaW52ZXJzZSBzaW5lICovXHJcblx0XHRjb25zdCBpID0gY29tcGxleCgwLCAxKTtcclxuXHRcdHJldHVybiB0aGlzLm11bHQoaSkuYWRkKGNvbXBsZXgoMSwgMCkuc3ViKHRoaXMuc3F1YXJlKCkpLnNxcnQoKSkubG4oKS5kaXYoaSk7XHJcblx0fVxyXG5cclxuXHRhdGFuKCkge1xyXG5cdFx0Y29uc3QgaSA9IGNvbXBsZXgoMCwgMSk7XHJcblx0XHRyZXR1cm4gQ29tcGxleC5kaXYoQ29tcGxleC5kaXYoaS5zdWIodGhpcyksIGkuYWRkKHRoaXMpKS5sbigpLCBpLnNjYWxlKDIpKTtcclxuXHR9XHJcblxyXG5cdHJvdGF0ZShhbmdsZSkge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhpcyBjb21wbGV4IG51bWJlciByb3RhdGVkIGJ5IGFuZ2xlIHJhZGlhbnMgKi9cclxuXHRcdHJldHVybiB0aGlzLm11bHQoKG5ldyBDb21wbGV4KDAsIGFuZ2xlKSkuZXhwKCkpO1xyXG5cdH1cclxuXHJcblx0ZG90KHopIHtcclxuXHRcdC8qIENvbXB1dGVzIHRoZSBFdWNsaWRlYW4gZG90IHByb2R1Y3Qgb2YgdGhlIGNvZWZmaWNpZW50cyBvZiB0aGlzIGNvbXBsZXggbnVtYmVyIGFuZCB6ICovXHJcblx0XHRyZXR1cm4gdGhpcy5yZSAqIHoucmUgKyB0aGlzLmltICogei5pbTtcclxuXHR9XHJcblxyXG5cdGFuZ2xlVG8oeikge1xyXG5cdFx0LyogQ29tcHV0ZXMgdGhlIGFuZ2xlIGJldHdlZW4gdGhpcyBjb21wbGV4IG51bWJlciBhbmQgeiAqL1xyXG5cdFx0LypcclxuXHRcdGFjb3MgdSp2L3V2ID0gdXZjb3ModClcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gTWF0aC5hY29zKHRoaXMuZG90KHopIC8gKHRoaXMubm9ybSgpICogei5ub3JtKCkpKTtcclxuXHR9XHJcblxyXG5cdHRvU3RyaW5nKCkge1xyXG5cdFx0LyogUmV0dXJucyB0aGUgc3RyaW5nIHJlcHJlc2VudGF0aW9uIG9mIHRoZSBjb21wbGV4IG51bWJlciBhcyBhbiBvcmRlcmVkIHBhaXIgKHJlKHopLCBpbSh6KSkgKi9cclxuXHRcdHJldHVybiBgKCR7dGhpcy5yZX0sJHt0aGlzLmltfSlgO1xyXG5cdH1cclxuXHJcblx0dG9MYXRleCgpIHtcclxuXHRcdC8qIFJldHVybnMgbGF0ZXggcmVwcmVzZW50YXRpb24gb2YgdGhlIGNvbXBsZXggbnVtYmVyICovXHJcblx0XHRjb25zdCByeCA9IHJvdW5kVG8odGhpcy5yZSwgMyksIHJ5ID0gcm91bmRUbyh0aGlzLmltLCAzKTtcclxuXHRcdHJldHVybiBgXFxcXGxlZnQoJHtyeH0sJHtyeX1cXFxccmlnaHQpYFxyXG5cdH1cclxuXHJcblx0ZXF1YWxzKHopIHtcclxuXHRcdC8qIFJldHVybnMgdHJ1ZSBpZmYgeiBlcXVhbHMgdGhpcyBjb21wbGV4IG51bWJlciwgZXhhY3RseSAqL1xyXG5cdFx0cmV0dXJuICh0aGlzLnJlID09IHoucmUgJiYgdGhpcy5pbSA9PSB6LmltKTtcclxuXHR9XHJcblxyXG5cdGVxdWFsc0Vwcyh6KSB7XHJcblx0XHQvKlxyXG5cdFx0UmV0dXJucyB0cnVlIGlmZiB6IGVxdWFscyB0aGlzIGNvbXBsZXggbnVtYmVyLCB3aXRoaW4gbnVtZXJpY2FsIHRvbGVyYW5jZSBFUFNJTE9OXHJcblx0XHRGb3IgZmxvYXRpbmcgcG9pbnQgcm91bmRpbmcgcHVycG9zZXNcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gKE1hdGguYWJzKHRoaXMucmUgLSB6LnJlKSA8IEVQU0lMT04gJiYgTWF0aC5hYnModGhpcy5pbSAtIHouaW0pIDwgRVBTSUxPTik7XHJcblx0fVxyXG5cclxuXHRtb2JpdXMoYSwgYiwgYywgZCkge1xyXG5cdFx0LypcclxuXHRcdEFwcGx5IHRoZSBNw7ZiaXVzIHRyYW5zZm9ybWF0aW9uIChheitiKS8oY3orZClcclxuXHRcdCovXHJcblx0XHRpZiAoQ29tcGxleC5pbmZpbml0ZSh0aGlzKSkge1xyXG5cdFx0XHRpZiAoYyAhPT0gMCkge1xyXG5cdFx0XHRcdHJldHVybiBDb21wbGV4LmRpdihhLCBjKTtcclxuXHRcdFx0fSBlbHNlIHtcclxuXHRcdFx0XHRyZXR1cm4gbmV3IENvbXBsZXgoSW5maW5pdHksIEluZmluaXR5KTtcclxuXHRcdFx0fVxyXG5cdFx0fSBlbHNlIHtcclxuXHRcdFx0Y29uc3QgZGVub21pbmF0b3IgPSBDb21wbGV4LmFkZChDb21wbGV4Lm11bHQoYywgdGhpcyksIGQpO1xyXG5cdFx0XHRpZiAoZGVub21pbmF0b3IgPT09IDApIHtcclxuXHRcdFx0XHRyZXR1cm4gbmV3IENvbXBsZXgoSW5maW5pdHksIEluZmluaXR5KTtcclxuXHRcdFx0fSBlbHNlIHtcclxuXHRcdFx0XHRjb25zdCBudW1lcmF0b3IgPSBDb21wbGV4LmFkZChDb21wbGV4Lm11bHQoYSwgdGhpcyksIGIpO1xyXG5cdFx0XHRcdHJldHVybiBDb21wbGV4LmRpdihudW1lcmF0b3IsIGRlbm9taW5hdG9yKTtcclxuXHRcdFx0fVxyXG5cdFx0fVxyXG5cdH1cclxuXHJcblx0c2luKCkge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgc2luZSBvZiB0aGUgY29tcGxleCBudW1iZXJcclxuXHRcdCovXHJcblx0XHRjb25zdCBpID0gY29tcGxleCgwLCAxKTtcclxuXHRcdGNvbnN0IHJvdGF0ZWQgPSB0aGlzLm11bHQoaSk7XHJcblx0XHRyZXR1cm4gcm90YXRlZC5leHAoKS5zdWIoIHJvdGF0ZWQuc2NhbGUoLTEpLmV4cCgpICkuZGl2KGkuc2NhbGUoMikpO1xyXG5cdH1cclxuXHJcblx0Y29zKCkge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgY29zaW5lIG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0Ki9cclxuXHRcdGNvbnN0IGkgPSBjb21wbGV4KDAsIDEpO1xyXG5cdFx0Y29uc3Qgcm90YXRlZCA9IHRoaXMubXVsdChpKTtcclxuXHRcdHJldHVybiByb3RhdGVkLmV4cCgpLmFkZCggcm90YXRlZC5zY2FsZSgtMSkuZXhwKCkgKS5zY2FsZSgwLjUpO1xyXG5cdH1cclxuXHJcblx0dGFuKCkge1xyXG5cdFx0LypcclxuXHRcdENhbGN1bGF0ZSB0aGUgdGFuZ2VudCBvZiB0aGUgY29tcGxleCBudW1iZXJcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gdGhpcy5zaW4oKS5kaXYodGhpcy5jb3MoKSlcclxuXHR9XHJcblxyXG5cdHNpbmgoKSB7XHJcblx0XHQvKlxyXG5cdFx0Q2FsY3VsYXRlIHRoZSBoeXBlcmJvbGljIHNpbmUgb2YgdGhlIGNvbXBsZXggbnVtYmVyXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHRoaXMuZXhwKCkuc3ViKHRoaXMuc2NhbGUoLTEpLmV4cCgpKS5zY2FsZSgwLjUpO1xyXG5cdH1cclxuXHJcblx0Y29zaCgpIHtcclxuXHRcdC8qXHJcblx0XHRDYWxjdWxhdGUgdGhlIGh5cGVyYm9saWMgY29zaW5lIG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB0aGlzLmV4cCgpLmFkZCh0aGlzLnNjYWxlKC0xKS5leHAoKSkuc2NhbGUoMC41KTtcclxuXHR9XHJcblxyXG5cdHRhbmgoKSB7XHJcblx0XHQvKlxyXG5cdFx0Q2FsY3VsYXRlIHRoZSBoeXBlcmJvbGljIHRhbmdlbnQgb2YgdGhlIGNvbXBsZXggbnVtYmVyXHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHRoaXMuc2luaCgpLmRpdih0aGlzLmNvc2goKSlcclxuXHR9XHJcblxyXG5cdHBvdyh6KSB7XHJcblx0XHQvKlxyXG5cdFx0Q2FsY3VsYXRlIHRoZSB6dGggcG93ZXIgb2YgdGhlIGNvbXBsZXggbnVtYmVyXHJcblx0XHQqL1xyXG5cdFx0aWYgKHRoaXMuZXF1YWxzRXBzKGNvbXBsZXgoMCwgMCkpKSB7XHJcblx0XHRcdHJldHVybiAoei5lcXVhbHNFcHMoY29tcGxleCgxLCAwKSkpID8gY29tcGxleCgxLCAwKSA6IGNvbXBsZXgoMCwgMCk7XHJcblx0XHR9XHJcblxyXG5cdFx0Y29uc3Qgc3ViQW5nID0gTWF0aC5hdGFuMih0aGlzLmltLCB0aGlzLnJlKTtcclxuXHRcdGNvbnN0IG5vcm1TcSA9IHRoaXMubm9ybVNxKCk7XHRcdFxyXG5cdFx0Y29uc3QgYW5nID0gMC41ICogei5pbSAqIE1hdGgubG9nKG5vcm1TcSkgKyB6LnJlICogc3ViQW5nO1xyXG5cdFx0Y29uc3Qgbm9ybSA9IE1hdGguZXhwKC16LmltICogc3ViQW5nKSAqIE1hdGgucG93KG5vcm1TcSwgMC41ICogei5yZSk7XHJcblx0XHRyZXR1cm4gY29tcGxleChub3JtICogTWF0aC5jb3MoYW5nKSwgbm9ybSAqIE1hdGguc2luKGFuZykpO1xyXG5cdH1cclxuXHJcblx0Y2xhbXAobWluLCBtYXgpIHtcclxuXHRcdC8qKlxyXG5cdFx0ICogY2xhbXAgdGhlIGNvbXBsZXggbnVtYmVyJ3Mgbm9ybSBiZXR3ZWVuIHR3byB2YWx1ZXNcclxuXHRcdCAqL1xyXG5cdFx0Y29uc3Qgbm9ybSA9IHoubm9ybSgpO1xyXG5cdFx0aWYgKCEobWluIDw9IG5vcm0gJiYgbm9ybSA8PSBtYXgpKSB7XHJcblx0XHRcdHJldHVybiB0aGlzLnNjYWxlKGNsYW1wKG5vcm0sIG1pbiwgbWF4KSAvIG5vcm0pO1xyXG5cdFx0fVxyXG5cdFx0cmV0dXJuIHRoaXM7XHJcblx0fVxyXG5cclxuXHRmcmFjKCkge1xyXG5cdFx0LyoqXHJcblx0XHQgKiByZXR1cm4gdGhlIGZyYWN0aW9uYWwgcGFydCBvZiBlYWNoIG9mIHRoZSBjb21wb25lbnRzIG9mIHRoZSBjb21wbGV4IG51bWJlclxyXG5cdFx0ICovXHJcblx0XHRyZXR1cm4gY29tcGxleChcclxuXHRcdFx0ZnJhYyh0aGlzLnJlKSxcclxuXHRcdFx0ZnJhYyh0aGlzLmltKSxcclxuXHRcdCk7XHJcblx0fVxyXG5cclxuXHQvKiAtLS0tLS0tLS0tIFN0YXRpYyBmdW5jdGlvbnMgLS0tLS0tLS0tLS0tLS0tLS0tLS0gKi9cclxuXHJcblx0c3RhdGljIG5vcm0oeikge1xyXG5cdFx0cmV0dXJuIHoubm9ybSgpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIG5vcm1TcSh6KSB7XHJcblx0XHRyZXR1cm4gei5ub3JtU3EoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBpbnYoeikge1xyXG5cdFx0cmV0dXJuIHouY29uaigpLnNjYWxlKDEgLyB6Lm5vcm1TcSgpKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBzY2FsZSh6LCBrKSB7XHJcblx0XHRyZXR1cm4gei5zY2FsZShrKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBhZGQoejEsIHoyKSB7XHJcblx0XHRyZXR1cm4gejEuYWRkKHoyKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBzdWIoejEsIHoyKSB7XHJcblx0XHRyZXR1cm4gejEuc3ViKHoyKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBtdWx0KHoxLCB6Mikge1xyXG5cdFx0cmV0dXJuIHoxLm11bHQoejIpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGRpdih6MSwgejIpIHtcclxuXHRcdHJldHVybiB6MS5kaXYoejIpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIHBvdyh6MSwgejIpIHtcclxuXHRcdHJldHVybiB6MS5wb3coejIpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGV4cCh6KSB7XHJcblx0XHRyZXR1cm4gei5leHAoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBsbih6KSB7XHJcblx0XHRyZXR1cm4gei5sbigpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIHNxcnQoeikge1xyXG5cdFx0cmV0dXJuIHouc3FydCgpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIHNpbih6KSB7XHJcblx0XHQvKlxyXG5cdFx0UmV0dXJuIHRoZSBzaW5lIG9mIHpcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gei5zaW4oKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBjb3Moeikge1xyXG5cdFx0LypcclxuXHRcdFJldHVybiB0aGUgY29zaW5lIG9mIHpcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gei5jb3MoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyB0YW4oeikge1xyXG5cdFx0LypcclxuXHRcdFJldHVybiB0aGUgdGFuZ2VudCBvZiB6XHRcdFxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB6LnRhbigpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGFzaW4oeikge1xyXG5cdFx0cmV0dXJuIHouYXNpbigpO1xyXG5cdH0gXHJcblxyXG5cdHN0YXRpYyBhY29zKHopIHtcclxuXHRcdHJldHVybiB6LmFjb3MoKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBhdGFuKHopIHtcclxuXHRcdHJldHVybiB6LmF0YW4oKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBzaW5oKHopIHtcclxuXHRcdC8qXHJcblx0XHRSZXR1cm4gdGhlIGh5cGVyYm9saWMgc2luZSBvZiB6XHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHouc2luaCgpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGNvc2goeikge1xyXG5cdFx0LypcclxuXHRcdFJldHVybiB0aGUgaHlwZXJib2xpYyBjb3NpbmUgb2YgelxyXG5cdFx0Ki9cclxuXHRcdHJldHVybiB6LmNvc2goKTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyB0YW5oKHopIHtcclxuXHRcdC8qXHJcblx0XHRSZXR1cm4gdGhlIGh5cGVyYm9saWMgdGFuZ2VudCBvZiB6XHJcblx0XHQqL1xyXG5cdFx0cmV0dXJuIHoudGFuaCgpO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGluZmluaXRlKHopIHtcclxuXHRcdC8qKlxyXG5cdFx0ICogcmV0dXJuIHdoZXRoZXIgb25lIG9yIGJvdGggb2YgdGhlIGNvbXBvbmVudHMgb2YgeiBpcyBpbmZpbml0ZSwgb3IgaWYgdGhlIG5vcm0gaXMgaW5maW5pdGVcclxuXHRcdCAqL1xyXG5cdFx0Y29uc3Qgbm9ybSA9IHoubm9ybSgpO1xyXG5cdFx0cmV0dXJuIChcclxuXHRcdFx0ei5yZSA9PT0gSW5maW5pdHkgfHwgei5yZSA9PT0gLUluZmluaXR5IFxyXG5cdFx0XHR8fCB6LmltID09PSBJbmZpbml0eSB8fCB6LmltID09PSAtSW5maW5pdHlcclxuXHRcdFx0fHwgbm9ybSA9PT0gSW5maW5pdHlcclxuXHRcdCk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgbmFuKHopIHtcclxuXHRcdC8qKlxyXG5cdFx0ICogcmV0dXJuIHdoZXRoZXIgb25lIG9yIGJvdGggY29tcG9uZW50cyBvZiB6IGlzIE5hTiBcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gaXNOYU4oei5yZSkgfHwgaXNOYU4oei5pbSk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgZ2FtbWEoeikge1xyXG5cdFx0aWYgKHoucmUgPCAwLjUpIHtcclxuXHRcdFx0Ly8gZ2FtbWEoMS16KWdhbW1hKHopID0gcGkgLyBzaW4ocGkqeilcclxuXHRcdFx0cmV0dXJuIENvbXBsZXguZGl2KCBjb21wbGV4KE1hdGguUEksIDAuMCksIENvbXBsZXgubXVsdCh6LnNjYWxlKE1hdGguUEkpLnNpbigpLCBDb21wbGV4LmdhbW1hKGNvbXBsZXgoMSAtIHoucmUsIC16LmltKSkpICk7XHJcblx0XHR9IGVsc2Uge1xyXG5cdFx0XHR6ID0gQ29tcGxleC5zdWIoeiwgY29tcGxleCgxLCAwKSk7IC8vIGFjY291bnQgZm9yIHN0dXBpZCBzaGlmdCBieSAxXHJcblx0XHRcdGxldCB4ID0gY29tcGxleChwVmFsdWVzWzBdLCAwKTtcclxuXHRcdFx0Zm9yIChsZXQgaT0xOyBpPHBWYWx1ZXMubGVuZ3RoOyBpKyspIHtcclxuXHRcdFx0XHR4ID0gQ29tcGxleC5hZGQoeCwgQ29tcGxleC5kaXYoIGNvbXBsZXgocFZhbHVlc1tpXSwgMC4wKSwgQ29tcGxleC5hZGQoeiwgY29tcGxleChpLCAwKSkgKSk7XHJcblx0XHRcdH1cclxuXHRcdFx0Y29uc3QgdCA9IENvbXBsZXguYWRkKHosIGNvbXBsZXgoNy41LCAwKSk7IC8vIGc9NywgZyswLjVcclxuXHRcdFx0cmV0dXJuIENvbXBsZXgucG93KHQsIHouYWRkKGNvbXBsZXgoMC41LCAwKSkpLm11bHQodC5zY2FsZSgtMSkuZXhwKCkpLm11bHQoeCkuc2NhbGUoTWF0aC5zcXJ0KDIgKiBNYXRoLlBJKSk7XHJcblx0XHR9XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgYmV0YSh6MSwgejIpIHtcclxuXHRcdHJldHVybiBDb21wbGV4LmRpdihDb21wbGV4Lm11bHQoQ29tcGxleC5nYW1tYSh6MSksIENvbXBsZXguZ2FtbWEoejIpKSwgQ29tcGxleC5nYW1tYSh6MS5hZGQoejIpKSk7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgbWF4KHoxLCB6Mikge1xyXG5cdFx0cmV0dXJuICh6MS5ub3JtU3EoKSA8IHoyLm5vcm1TcSgpKSA/IHoyIDogejE7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgbWluKHoxLCB6Mikge1xyXG5cdFx0cmV0dXJuICh6MS5ub3JtU3EoKSA+IHoyLm5vcm1TcSgpKSA/IHoyIDogejE7XHJcblx0fVxyXG5cclxuXHRzdGF0aWMgY2xhbXAoeiwgbWluLCBtYXgpIHtcclxuXHRcdHJldHVybiB6LmNsYW1wKG1pbiwgbWF4KTtcclxuXHR9XHJcblxyXG5cdHN0YXRpYyBmcmFjKCkge1xyXG5cdFx0cmV0dXJuIHouZnJhYygpO1xyXG5cdH1cclxuXHJcblx0LyogLS0tLS0tLS0tLS0tLS0tIEluLXBsYWNlIG9wZXJhdGlvbnMgLS0tLS0tLS0tLS0tLS0tLS0tLS0tICovXHJcblxyXG5cdGlhZGQoeikge1xyXG5cdFx0dGhpcy5yZSArPSB6LnJlO1xyXG5cdFx0dGhpcy5pbSArPSB6LmltO1xyXG5cdH1cclxuXHJcblx0aXN1Yih6KSB7XHJcblx0XHR0aGlzLnJlIC09IHoucmU7XHJcblx0XHR0aGlzLmltIC09IHouaW07XHJcblx0fVxyXG5cclxuXHJcbn1cclxuXHJcbmZ1bmN0aW9uIGNvbXBsZXgocmVhbCwgaW1hZ2luYXJ5KSB7XHJcblx0LyogaW5zdGFudGlhdGUgYSBDb21wbGV4IHdpdGhvdXQgbmV3IGtleXdvcmQgKi9cclxuXHRyZXR1cm4gbmV3IENvbXBsZXgocmVhbCwgaW1hZ2luYXJ5KTtcclxufVxyXG5cclxuXHJcbmZ1bmN0aW9uIGludGVncmF0ZU92ZXJQYXJhbWV0ZXIoZiwgYT0wLCBiPTEsIHNhbXBsZXM9MTAwKSB7XHJcblx0LypcclxuXHRyZXR1cm4gdGhlIGxpbmUgaW50ZWdyYWwgb2YgdGhlIENvbXBsZXgtdmFsdWVkIGZ1bmN0aW9uIGZcclxuXHRwYXJhbWV0ZXJpemVkIGJ5IHRoZSByZWFsIHBhcmFtZXRlciB0IG9uIFthLCBiXSxcclxuXHR1c2luZyB0aGUgc3BlY2lmaWVkIG51bWJlciBvZiBzYW1wbGVzXHJcblx0Ki9cclxuXHRjb25zdCBzdGVwID0gKGIgLSBhKSAvIHNhbXBsZXM7XHJcblx0bGV0IHQgPSBhICsgc3RlcCAvIDI7IC8vIG1pZHBvaW50IGFwcHJveGltYXRpb25cclxuXHRsZXQgcmVzdWx0ID0gY29tcGxleCgwLCAwKTtcclxuXHRmb3IgKGxldCBpPTA7IGk8c2FtcGxlczsgaSsrKSB7XHJcblx0XHRyZXN1bHQuaWFkZChmKHQpLnNjYWxlKHN0ZXApKTtcclxuXHRcdHQgKz0gc3RlcDtcclxuXHR9XHJcblx0cmV0dXJuIHJlc3VsdDtcclxufVxyXG5cclxuZnVuY3Rpb24gaW50ZWdyYXRlT3ZlckNvbXBsZXhWYXJpYWJsZShmLCB6MCwgejEsIHNhbXBsZXM9MTAwKSB7XHJcblx0LypcclxuXHRyZXR1cm4gdGhlIGxpbmUgaW50ZWdyYWwgb2YgdGhlIENvbXBsZXgtdmFsdWVkIGZ1bmN0aW9uIGZcclxuXHRvbiB0aGUgbGluZSBiZXR3ZWVuIHowIGFuZCB6MSwgdXNpbmcgdGhlIHNwZWNpZmllZCBudW1iZXIgb2YgaW50ZXJ2YWxzXHJcblx0Ki9cclxuXHRjb25zdCBzdGVwID0gQ29tcGxleC5zdWIoejEsIHowKS5zY2FsZSgxIC8gc2FtcGxlcyk7XHJcblx0Y29uc3Qgc3RlcE5vcm0gPSBzdGVwLm5vcm0oKTtcclxuXHRsZXQgeiA9IENvbXBsZXguYWRkKHowLCBzdGVwLnNjYWxlKDAuNSkpO1xyXG5cdGxldCByZXN1bHQgPSBjb21wbGV4KDAsIDApO1xyXG5cdGZvciAobGV0IGk9MDsgaTxzYW1wbGVzOyBpKyspIHtcclxuXHRcdHJlc3VsdC5pYWRkKGYoeikuc2NhbGUoc3RlcE5vcm0pKTtcclxuXHRcdHouaWFkZChzdGVwKTtcclxuXHR9XHJcblx0cmV0dXJuIHJlc3VsdDtcclxufVxyXG5cclxuXHJcblxyXG5jbGFzcyBWZWN0b3Ige1xyXG5cclxuXHQvKlxyXG5cdFJlcHJlc2VudHMgYSB2ZWN0b3IgaW4gQ15uXHJcblx0Ki9cclxuXHJcblx0Y29uc3RydWN0b3IodmFsdWVzKSB7XHJcblx0XHR0aGlzLnZhbHVlcyA9IHZhbHVlcztcclxuXHRcdHRoaXMuZGltZW5zaW9uID0gdmFsdWVzLmxlbmd0aDtcclxuXHR9XHJcblxyXG5cdGdldChpbmRleCkge1xyXG5cdFx0cmV0dXJuIHRoaXMudmFsdWVzW2luZGV4XTtcclxuXHR9XHJcblxyXG5cdHNjYWxlKHopIHtcclxuXHRcdC8qXHJcblx0XHRTY2FsZSBlYWNoIGNvbXBvbmVudCBieSB0aGUgY29tcGxleCBudW1iZXIgelxyXG5cdFx0Ki9cclxuXHRcdGNvbnN0IHJlc3VsdCA9IFtdO1xyXG5cdFx0Zm9yIChsZXQgaT0wOyBpPHRoaXMuZGltZW5zaW9uOyBpKyspIHtcclxuXHRcdFx0cmVzdWx0LnB1c2goQ29tcGxleC5tdWx0KHRoaXMuZ2V0KGkpLCB6KSk7XHJcblx0XHR9XHJcblx0XHRyZXR1cm4gdmVjdG9yKHJlc3VsdCk7XHJcblx0fVxyXG5cclxuXHRyZWFsU2NhbGUoaykge1xyXG5cdFx0LypcclxuXHRcdFNjYWxlIGVhY2ggY29tcG9uZW50IGJ5IHRoZSByZWFsIG51bWJlciBrXHJcblx0XHQqL1xyXG5cdFx0Y29uc3QgcmVzdWx0ID0gW107XHJcblx0XHRmb3IgKGxldCBpPTA7IGk8dGhpcy5kaW1lbnNpb247IGkrKykge1xyXG5cdFx0XHRyZXN1bHQucHVzaCh0aGlzLmdldChpKS5zY2FsZShrKSk7XHJcblx0XHR9XHJcblx0XHRyZXR1cm4gdmVjdG9yKHJlc3VsdCk7XHJcblx0fVxyXG5cclxuXHRhZGQoWikge1xyXG5cdFx0Y29uc3QgcmVzdWx0ID0gW107XHJcblx0XHRmb3IgKGxldCBpPTA7IGk8dGhpcy5kaW1lbnNpb247IGkrKykge1xyXG5cdFx0XHRyZXN1bHQucHVzaChDb21wbGV4LmFkZCh0aGlzLmdldChpKSwgWi5nZXQoaSkpKTtcclxuXHRcdH1cclxuXHRcdHJldHVybiB2ZWN0b3IocmVzdWx0KTtcclxuXHR9XHJcblx0XHJcblx0ZG90KFopIHtcclxuXHRcdC8qXHJcblx0XHRDb21wdXRlIChIZXJtaXRpYW4pIGRvdCBwcm9kdWN0XHJcblx0XHRFcXVpdmFsZW50IHRvIHN0YW5kYXJkIGRvdCBwcm9kdWN0IGluIHRoZSBjYXNlIG9mIHJlYWwgdmVjdG9yc1xyXG5cdFx0Ki9cclxuXHRcdGNvbnN0IHJlc3VsdCA9IGNvbXBsZXgoMCwgMCk7XHJcblx0XHRmb3IgKGxldCBpPTA7IGk8dGhpcy5kaW1lbnNpb247IGkrKykge1xyXG5cdFx0XHRyZXN1bHQuaWFkZChDb21wbGV4Lm11bHQodGhpcy5nZXQoaSksIFouZ2V0KGkpLmNvbmooKSkpO1xyXG5cdFx0fVxyXG5cdFx0cmV0dXJuIHJlc3VsdDtcclxuXHR9XHJcblxyXG5cdG5vcm0oKSB7XHJcblx0XHQvKlxyXG5cdFx0Q29tcHV0ZSBub3JtIG9mIHRoZSB2ZWN0b3JcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gTWF0aC5zcXJ0KHRoaXMuZG90KHRoaXMpLnJlKTtcclxuXHR9XHJcblxyXG5cdGFuZ2xlVG8oWikge1xyXG5cdFx0LypcclxuXHRcdENvbXB1dGUgdGhlIGFuZ2xlIGJldHdlZW4gdGhlIHZlY3RvciBhbmQgdmVjdG9yIFpcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gTWF0aC5hY29zKCB0aGlzLmRvdChaKS5yZSAvICh0aGlzLm5vcm0oKSAqIFoubm9ybSgpKSApO1xyXG5cdH1cclxuXHJcblx0c3RhdGljIGxlcnAoWjEsIFoyLCB0KSB7XHJcblx0XHQvKlxyXG5cdFx0TGluZWFybHkgaW50ZXJwb2xhdGVcclxuXHRcdCovXHJcblx0XHRyZXR1cm4gWjEucmVhbFNjYWxlKDEgLSB0KS5hZGQoWjIuc2NhbGUodCkpO1xyXG5cdH1cclxuXHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiB2ZWN0b3IoLi4udmFsdWVzKSB7XHJcblx0LypcclxuXHRVdGlsaXR5IGZ1bmN0aW9uIGZvciBlYXN5IGluc3RhbnRpYXRpb24gb2YgVmVjdG9yIG9iamVjdHMuXHJcblx0QWxsb3dzIGNvbXBvbmVudHMgdG8gYmUgcGFzc2VkIGluZGl2aWR1YWxseSBvciBpbiBhbiBhcnJheVxyXG5cdCovXHJcblx0aWYgKHZhbHVlc1swXSBpbnN0YW5jZW9mIENvbXBsZXgpIHtcclxuXHRcdHJldHVybiBuZXcgVmVjdG9yKHZhbHVlcyk7XHJcblx0fSBlbHNlIHtcclxuXHRcdHJldHVybiBuZXcgVmVjdG9yKHZhbHVlc1swXSk7XHJcblx0fVxyXG59XHJcblxyXG5cclxubW9kdWxlLmV4cG9ydHMgPSB7XHJcblx0Q29tcGxleCwgY29tcGxleCxcclxuXHRpbnRlZ3JhdGVPdmVyQ29tcGxleFZhcmlhYmxlLCBpbnRlZ3JhdGVPdmVyUGFyYW1ldGVyLFxyXG5cdFZlY3RvciwgdmVjdG9yLFxyXG59OyJdfQ==
