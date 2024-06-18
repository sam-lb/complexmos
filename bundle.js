(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
/**
 * I think per-pixel function evaluation in javascript is fast enough to run in real time; rendering is the bottleneck
 * so JS generates the mesh and coloring can be done with a dynamic property passed as a uniform in reGL
 */



const { valueScope } = require("./scope.js");
const { tracker } = require("./parsing/errors.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("./parsing/pratt/expressions.js");
const { TokenType } = require("./parsing/pratt/tokentype.js");
const { complex, Complex } = require("./math/complex.js");




function evaluate(ast) {
    // at this point, the ast has gone through the lexer and parser's checks already
    if (ast instanceof AssignExpression) {
        // the expression is trying to define a function
        const left = ast.mLeft;
        if (left instanceof CallExpression) {
            // the expression is defining a function
            valueScope[left.mFunction.mName] = new Evaluatable(ast.mRight,
                left.mArgs.map((arg) => (arg instanceof OperatorExpression) ? arg.mLeft.mName : arg.mName)
            );
        } else {
            // the expression is defining a variable
            valueScope[left.mName] = new Evaluatable(ast.mRight);            
        }
        return null;
    } else {
        return new Evaluatable(ast);
    }
}


class Evaluatable {

    constructor(ast, args=null) {
        this.ast = ast;
        this.args = (args === null) ? [] : args;
    }

    call(args) {
        for (const arg of this.args) {
            if (args[arg] === undefined) {
                tracker.error(`function requires argument ${arg}`);
                return null;
            }
        }
        try {
            return this._call(this.ast, args);
        } catch (error) {
            tracker.error(`An unrecognized error has occurred: ${error.message}`);
            return null;
        }
    }

    _call(ast, args=null) {
        args = (args === null) ? {} : args;
        if (ast instanceof OperatorExpression) {
            const arg1 = this._call(ast.mLeft, args);
            const arg2 = this._call(ast.mRight, args);

            switch(ast.mOperator) {
                case (TokenType.PLUS):
                    return Complex.add(arg1, arg2);        
                case (TokenType.MINUS):
                    return Complex.sub(arg1, arg2);
                case (TokenType.ASTERISK):
                    return Complex.mult(arg1, arg2);
                case (TokenType.SLASH):
                    return Complex.div(arg1, arg2);
                case (TokenType.CARET):
                    return Complex.pow(arg1, arg2);
                default:
                    tracker.error("unexpected operator encountered");
                    return null;
            }
        } else if (ast instanceof PrefixExpression) {
            const arg1 = this._call(ast.mRight, args);
            if (ast.mOperator === TokenType.MINUS) {
                return arg1.scale(-1);
            } else {
                tracker.error("unexpected prefix operator encountered");
            }
        } else if (ast instanceof CallExpression) {
            const functionArgs = ast.mArgs.map(arg => this._call(arg, args));
            // we don't check local argument scope here because higher order functions aren't allowed
            if (valueScope[ast.mFunction] !== undefined) {
                const func = valueScope[ast.mFunction];
                if (func instanceof Evaluatable) {
                    const argMap = {};
                    func.args.forEach((key, index) => argMap[key] = functionArgs[index]);
                    return func.call(argMap);
                } else {
                    return func(...functionArgs);
                }
            } else {
                tracker.error(`could not resolve function ${ast.mFunction}. Note: higher order functions are not yet supported`);
            }
        } else if (ast instanceof NameExpression) {
            if (args[ast.mName] !== undefined) {
                return args[ast.mName];
            } else if (valueScope[ast.mName] !== undefined) {
                let value = valueScope[ast.mName];
                if (value instanceof Evaluatable) value = value.call();
                return value;
            } else {
                tracker.error(`could not resolve variable ${ast.mName}`);
                return null;
            }
        } else if (ast instanceof NumberExpression) {
            return complex(ast.mNumber, 0);
        } else {
            tracker.error("what is even going on");
        }
    }
}



module.exports = {
    evaluate, Evaluatable,
};
},{"./math/complex.js":2,"./parsing/errors.js":9,"./parsing/pratt/expressions.js":12,"./parsing/pratt/tokentype.js":18,"./scope.js":20}],2:[function(require,module,exports){
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
		this.add(this.square().sub(new Complex(1, 0)).sqrt()).ln().div(complex(0, 1));
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
},{}],3:[function(require,module,exports){

const { Euclid } = require("./geometry.js");
const { complex, Complex, integrateOverParameter } = require("./complex.js");




function parameterizePoints(points) {
    /*
    Return a parameteric function consisting of straight lines connecting points
    Note: does NOT necessarily produce a constant time parameterization
    */
    const N = points.length;
    const segmentLength = 1 / (N - 1);
    return (t) => {
        if (t === 1) {
            return points[N - 1];
        } else {
            const index = Math.floor(t / segmentLength);
            return Euclid.lerp(points[index], points[index+1], t / segmentLength - Math.floor(t / segmentLength));
        }
    };
}



function fourierCoefficients(f, N) {
    /*
    Calculate the fourier coefficients for the complex function f 
    parameterized on [0, 1]. 
    Coefficients from -N, -N+1, ... , N-1, N
    */
    const coefs = [];
    for (let n=-N; n<=N; n++) {
        const F = (t) => {
            return Complex.mult(f(t), Complex.exp( complex(0, -2 * Math.PI * n * t) ));
        };
        coefs.push(
            integrateOverParameter(F, 0, 1)
        );
    }
    return coefs;
}


function fourierSeries(f, N) {
    const coefs = fourierCoefficients(f, N);
    return (t) => {
        const result = complex(0, 0);
        for (let n=-N; n<=N; n++) {
            const coef = coefs[n + N];
            result.iadd(
                Complex.mult(coef, Complex.exp( complex(0, 2 * Math.PI * n * t) ))
            );
        }
        return result;
    };
}


module.exports = {
    parameterizePoints,
    fourierCoefficients,
    fourierSeries,
};
},{"./complex.js":2,"./geometry.js":4}],4:[function(require,module,exports){

const { Complex, complex } = require("./complex.js");


class Euclid {

	/*
	Collection of Euclidean geometry functions
	*/

	static lineIntersection(p1, v1, p2, v2) {
		/*
		Computes the intersection of the lines p1 + v1 * t and p2 + v2 * t.
		Returns null if the lines do not intersect.
		*/
		const x1 = p1.re, y1 = p1.im;
		const x2 = p1.re + v1.re, y2 = p1.im + v1.im;
		const x3 = p2.re, y3 = p2.im;
		const x4 = p2.re + v2.re, y4 = p2.im + v2.im;

		const denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
		if (denom == 0) return null;

		const t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
		return p1.add(v1.scale(t));
	}

	static midpoint(p1, p2) {
		/*
		Computes the midpoint of p1 and p2
		*/
		return p1.add(p2).scale(0.5);
	}

	static lerp(p1, p2, t) {
		/*
		Computes linear interpolation between p1 and p2
		t is on [0, 1]: lerp(p1, p2, 0) = p1, lerp(p1, p2, 1) = p2
		*/
		return Complex.add( p1.scale(1 - t), p2.scale(t) );
	}

	static circleCenter(p1, p2, p3) {
		/*
		Computes the center of the circle passing through the three points p1, p2, p3
		*/
		return Euclid.lineIntersection(Euclid.midpoint(p1, p2), p2.sub(p1).perp(),
										Euclid.midpoint(p2, p3), p3.sub(p2).perp());
	}

	static centroid(P) {
		/*
		Computes the centroid of the points in P
		*/
		let xTotal = 0, yTotal = 0;
		for (let point of P) {
			xTotal += point.re;
			yTotal += point.im;
		}
		return complex(xTotal / P.length, yTotal / P.length);
	}

	static distance(p1, p2) {
		/*
		Computes the distance between p1 and p2
		*/
		return p2.sub(p1).norm();
	}

	static project(p1, p2) {
		return p2.scale(p1.dot(p2) / p2.normSq());
	}

}


class Poincare {

	/* Collection of functions for computations in the poincare disk model of the hyperbolic plane */

	static translatePToOrigin(z, P) {
		/*
		Computes a mobius transformation on z that takes P to 0 and preserves the unit disk
		*/
		return z.sub(P).div(complex(1, 0).sub(P.conj().mult(z)));
	}

	static translateOriginToP(z, P) {
		/*
		Computes a mobius transformation on z taking 0 to P and preserves the unit disk
		(inverse of translatePToOrigin) */
		return z.add(P).div(complex(1, 0).add(P.conj().mult(z)));
	}

	static segment(t, A, B) {
		/*
		Evaluates a parameterization of the geodesic segment between A and B at time t.
		segment(0, A, B) = A.
		segment(1, A, B) = B.
		*/
		return Poincare.translateOriginToP(Poincare.translatePToOrigin(B, A).scale(t), A);
	}

	static line(t, A, B) {
		/*
		Evaluates a parameterization of the geodesic through A and B at time t.
		line(0, A, B) = start of line (beginning point on circle at infinity)
		line(1, A, B) = end of line (terminal point on circle at infinity)
		Note that line(0, A, B) and line(1, A, B) are not actually points on the geodesic.
		*/
		return Poincare.translateOriginToP(Poincare.translatePToOrigin(B, A).unit().scale(2 * t - 1), A);	
	}

	static regPolyDist(p, q) {
		/*
		Computes the (Euclidean) distance to vertices of a regular p-gon with interior
		angle 2*pi/q (for (p, q) tessellation).
		Note: (p-2) * (q-2) must be greater than 4
		*/
		if ((p-2) * (q-2) <= 4) {
			console.error(`Error: cannot compute regular polygon distance for p=${p}, q=${q}`);
			return;
		}

		const tan1 = Math.tan(Math.PI / 2 - Math.PI / q);
		const tan2 = Math.tan(Math.PI / p);
		return Math.sqrt((tan1 - tan2) / (tan1 + tan2));
	}

	static polygon(N, verts) {
		// make a hyperbolic polygon with as close to N points as possible while guaranteeing all corners
		if (verts.length < 2) {
			console.error("Error: can't draw polygon with less than 2 points")
		}
		const result = [];
		const nPerSide = Math.ceil(N / verts.length);
		
		verts = verts.slice();
		verts.push(verts[0]);
		for (let i=0; i<verts.length-1; i++) {
			const space = linspace(0, 1, nPerSide);
			for (let value of space) {
				result.push(Poincare.segment(value, verts[i], verts[i+1]));
			}
		}
		return result;
	}

	static rotate(z, P, angle) {
		/* Computes the hyperbolic rotation of z about P by angle radians */
		return Poincare.translateOriginToP(Poincare.translatePToOrigin(z, P, true).rotate(angle), P, true);
	}

	static rotateMultiple(Z, P, angle) {
		/* Computes the hyperbolic rotation of all the points in Z about P by angle radians */
		const result = [];
		for (let vert of Z) {
			result.push(this.rotate(vert, P, angle));
		}
		return result;
	}

	static unitCircleInvert(z) {
		/* Computes the inversion of z through the unit circle */
		return z.conj().inv();
	}

	static circleInvert(z, r, P) {
		/* Computes the inversion of z through the circle of radius r centered at P */
		return P.add(complex(r * r, 0).div(z.sub(P).conj()));
	}

	static reflect(z, p1, p2) {
		/* Computes the inversion of z through the geodesic passing through p1 and p2 */
		const center = Euclid.circleCenter(p1, p2, Poincare.unitCircleInvert(p1));

		if (center == null || center.norm() > 1000 || isNaN(center.re) || isNaN(center.im)) {
			// the points are presumably on a radial line through the origin
			return p1.add(Euclid.project(z.sub(p1), p2.sub(p1))).scale(2).sub(z);
		}
		return Poincare.circleInvert(z, Euclid.distance(p1, center), center);
	}

	static reflectMultiple(Z, p1, p2) {
		/* Computes the inversion of all points in Z through the geodesic passing through p1 and p2 */
		const result = [];
		for (let vert of Z) {
			result.push(this.reflect(vert, p1, p2));
		}
		return result;
	}

	static inverseCayley(z) {
		/*
		Inverse of the Cayley transform (map from upper half plane to unit disk)
		See https://en.wikipedia.org/wiki/Cayley_transform#Complex_homography
		and https://www.desmos.com/calculator/ucug6yw6bh
		*/
		const one = complex(1, 0);
		return one.add(z).div(one.sub(z)).mult(complex(0, 1));
	}

	static toKleinDisk(z) {
		/*
		Take a point z in the Poincaré disk to the Klein disk
		(Inverse stereographic projection to hemisphere then orthogonal projection to disk)
		*/
		const denom = 0.5 * (1 + z.normSq());
		return complex(z.re / denom, z.im / denom);
	}

	static hypDistance(z1, z2) {
		/*
		Given two points, compute the hyperbolic distance between them according to the Poincare metric
		*/
		const z = z1.sub(z2);
		const euclideanDistanceToOrigin = z.norm();
		return Math.log((1 + euclideanDistanceToOrigin) / (1 - euclideanDistanceToOrigin));
	}

}


module.exports = {
	Euclid, Poincare,
};
},{"./complex.js":2}],5:[function(require,module,exports){




function snormal(X) {
    const norm = Math.sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
    return [
        X[0] / norm,
        X[1] / norm,
        X[2] / norm,
    ];
}

function ssub(s, X, Y) {
    return [
        X[0] + s * Y[0],
        X[1] + s * Y[1],
        X[2] + s * Y[2],
    ];
}

function sscale(s, X) {
    return [
        s * X[0],
        s * X[1],
        s * X[2],
    ];
}



function icosphere(subdivisions=0) {

    // original icosahedron data

    const scale = 1 / (2 * Math.sin(2 * Math.PI / 5));
    const phi = scale * (1 + Math.sqrt(5)) / 2;

    const vertices = [
        [-scale, phi, 0], [scale, phi, 0], [-scale, -phi, 0], [scale, -phi, 0],
        [0, -scale, phi], [0, scale, phi], [0, -scale, -phi], [0, scale, -phi],
        [phi, 0, -scale], [phi, 0, scale], [-phi, 0, -scale], [-phi, 0, scale],
    ];

    const triangleIndices = [
        [0, 11, 5],[0, 5, 1],[0, 1, 7],[0, 7, 10],[0, 10, 11],
        [11, 10, 2],[5, 11, 4],[1, 5, 9],[7, 1, 8],[10, 7, 6],
        [3, 9, 4],[3, 4, 2],[3, 2, 6],[3, 6, 8],[3, 8, 9],
        [9, 8, 1],[4, 9, 5],[2, 4, 11],[6, 2, 10],[8, 6, 7],
    ];

    let triangles = [];
    for (let index of triangleIndices) {
        triangles.push([vertices[index[0]], vertices[index[1]], vertices[index[2]]]);
    }

    // split into 9:

    // const baryCoefs = [
    //     [0, 0], [1/3, 0], [2/3, 0], [1, 0], [0, 1/3],
    //     [1/3, 1/3], [2/3, 1/3], [0, 2/3], [1/3, 2/3], [0, 1],
    // ];

    // const trindices = [
    //     [0, 1, 4], [1, 2, 5], [2, 3, 6], [4, 5, 1], [5, 6, 2],
    //     [4, 5, 7], [5, 6, 8], [7, 8, 5], [7, 8, 9],
    // ];

    // split into 4:

    const baryCoefs = [
        [0, 0], [.5, 0], [1, 0],
        [0, .5], [.5, .5], [0, 1],
    ];

    const trindices = [
        [0, 1, 3], [1, 2, 4],
        [1, 3, 4], [3, 4, 5],
    ];

    for (let i=0; i<subdivisions; i++) {
        const newTriangles = [];
        for (let tri of triangles) {
            let newVertices = [];
            const u = [tri[1][0] - tri[0][0], tri[1][1] - tri[0][1], tri[1][2] - tri[0][2]];
            const v = [tri[2][0] - tri[0][0], tri[2][1] - tri[0][1], tri[2][2] - tri[0][2]];

            for (let coef of baryCoefs) {
                newVertices.push(ssub(coef[1], ssub(coef[0], tri[0], u), v));
            }

            for (let trindex of trindices) {
                newTriangles.push([
                    newVertices[trindex[0]], newVertices[trindex[1]], newVertices[trindex[2]]
                ]);
            }
        }
        triangles = newTriangles.slice();
    }
    
    // normalize vertices
    for (let i=0; i<triangles.length; i++) {
        triangles[i][0] = snormal(triangles[i][0]);
        triangles[i][1] = snormal(triangles[i][1]);
        triangles[i][2] = snormal(triangles[i][2]);
    }
    return triangles;

}


module.exports = {
    snormal, ssub, sscale, // remove these later (Fix this hacky garbage)
    icosphere,
};
},{}],6:[function(require,module,exports){
const config = {
	EPSILON: 1e-8, // acceptable numerical error (shivers)
};

const cheapNumberRound = (x) => {
	// check if x is within EPSILON of floor(x) or ceil(x), if so return the one of these it's close to
	let u = Math.ceil(x), l = Math.floor(x);
	if (Math.abs(x-u) < config.EPSILON) return u;
	if (Math.abs(x-l) < config.EPSILON) return l;
	return x;
};


class Matrix {

	constructor(arr) {
		if (typeof arr[0] === "number") {
			arr = [arr];
		}
		this.mat = arr;
		this.cols = this.mat[0].length;
		this.rows = this.mat.length;
		if (this.rows === 0 || this.cols === 0) throw new Error(`Invalid matrix dimensions ${this.rows} x ${this.cols}.`);
		this.isSquare = this.rows === this.cols;
	}

	static randomMatrix(m, n) {
		// generate a random m by n matrix
		let arr = [];
		for (let i=0; i<m; i++) {
			arr[i] = [];
			for (let j=0; j<n; j++) {
				arr[i][j] = Math.random();
			}
		}
		return new Matrix(arr);
	}

	static randomIntegerMatrix(m, n, lower=0, upper=100) {
		// generate a random m by n matrix of integers between lower and upper
		lower = Math.floor(lower);
		upper = Math.floor(upper+1);
		let arr = [], range = upper - lower;
		for (let i=0; i<m; i++) {
			arr[i] = [];
			for (let j=0; j<n; j++) {
				arr[i][j] = lower + Math.floor(Math.random() * range);
			}
		}
		return new Matrix(arr);
	}

	static emptyMatrix(m, n) {
		// generate an m by n matrix of zeros
		return Matrix.fillMatrix(m, n, 0);
	}

	static fillMatrix(m, n, v) {
		// generate an m by n matrix, all entries having value v
		const arr = [];
		for (let i=0; i<m; i++) {
			arr[i] = [];
			for (let j=0; j<n; j++) {
				arr[i][j] = v;
			}
		}
		return new Matrix(arr);
	}

	static identity(n) {
		// generate the n by n identity matrix
		const arr = [];
		for (let i=0; i<n; i++) {
			arr[i] = [];
			for (let j=0; j<n; j++) {
				if (i === j) {
					arr[i][j] = 1;
				} else {
					arr[i][j] = 0;
				}
			}
		}
		return new Matrix(arr);
	}

	static determinant2x2(a, b, c, d) {
		/* calculate the determinant of the matrix:
			| a   b |
			| c   d |
		This is a separate function because it's frequently used
		*/
		return a * d - b * c;
	}

	static rotationMatrix2D(angle) {
		/** Generate a 2d rotation matrix for the given angle */
		return new Matrix([[Math.cos(angle), -Math.sin(angle)], [Math.sin(angle), Math.cos(angle)]]);
	}

	static rotationMatrix3D(angleX, angleY, angleZ) {
		/** Generate a 3d rotation matrix for the given Euler angles */
		const sinX = Math.sin(angleX), cosX = Math.cos(angleX);
		const sinY = Math.sin(angleY), cosY = Math.cos(angleY);
		const sinZ = Math.sin(angleZ), cosZ = Math.cos(angleZ);
		const rotX = new Matrix([
			[1, 0, 0],
			[0, cosX, -sinX],
			[0, sinX, cosX],
		]);
		const rotY = new Matrix([
			[cosY, 0, sinY],
			[0, 1, 0],
			[-sinY, 0, cosY],
		]);
		const rotZ = new Matrix([
			[cosZ, -sinZ, 0],
			[sinZ, cosZ, 0],
			[0, 0, 1],
		]);

		return Matrix.multiply(Matrix.multiply(rotX, rotY), rotZ);
	}

	static multiply(mat1, mat2) {
		if (mat1.cols !== mat2.rows) throw new Error("Incompatible matrix dimensions. must be m x r and r x n");
		const res = [];
		for (let i=0; i<mat1.rows; i++) {
			const row = [];
			for (let j=0; j<mat1.cols; j++) {
				let dot = 0;
				for (let k=0; k<mat2.rows; k++) {
					dot += mat1.mat[i][k] * mat2.mat[k][j];
				}
				row.push(dot);
			}
			res.push(row);
		}
		return new Matrix(res);
	}

	static scale(mat, k) {
		let newMat = mat.copy();
		newMat.scale(k);
		return newMat;
	}

	static add(mat1, mat2) {
		let newMat = mat1.copy();
		newMat.add(mat2);
		return newMat;
	}

	static sub(mat1, mat2) {
		let newMat = mat1.copy();
		newMat.sub(mat2);
		return newMat;
	}

	static detCofactorExpansion(mat) {
		// compute the determinant of mat using cofactor expansion
		if (mat.isSquare) {
			if (mat.rows > 10) { throw new Error("For larger matrices, use Matrix.det (cofactor expansion is inefficient and computation time grows with the factorial of matrix dimension)."); }
			if (mat.rows === 2 && mat.cols === 2) {
				return Matrix.determinant2x2(mat.mat[0][0], mat.mat[0][1], mat.mat[1][0], mat.mat[1][1]);
			} else {
				let sign = 1, bottom = mat.withoutRow(0), res = 0;
				for (let j=0; j<mat.rows; j++) {
					res += sign * mat.mat[0][j] * Matrix.detCofactorExpansion(bottom.withoutColumn(j));
					sign = -sign;
				}
				return res;
			}
		}
		throw new Error("The determinant is only defined for square matrices.");
	}

	static det(mat) {
		// compute the determinant of mat, more efficiently.

		/*
			For a triangular matrix, the determinant is the product of the entries along the diagonal.
			Row-reduce to get into upper triangular form, then compute the product.
		*/

		let prod = null;
		if (mat.isSquare) {
			let newMat = Matrix.LUDecomposition(mat).U;
			for (let d=0; d<newMat.rows; d++) {
				if (prod === null) {
					prod = newMat.mat[d][d];
				} else {
					prod *= newMat.mat[d][d];
				}
				if (prod === 0) break;
			}
			return prod;
		}
		throw new Error("The determinant is only defined for square matrices.");
	}

	static ref(mat) {
		// convert mat to row-echelon form
		let newMat = mat.copy();
		newMat.ref();
		return newMat;
	}

	static rref(mat) {
		// return mat reduced row-echelon form
		let newMat = mat.copy();
		newMat.rref();
		return Matrix.cheapMatrixRound(newMat);
	}

	static columnNorm(mat, j) {
		let squaredTotal = 0;
		for (let i=0; i<mat.rows; i++) {
			squaredTotal += mat.get(i, j) * mat.get(i, j);
		}
		return Math.sqrt(squaredTotal);
	}

	static householderReflection(mat) {
		/**
		 * compute the householder reflection matrix for v
		 * v should be of shape (n, 1)
		 */
		if (mat.cols !== 1) throw new Error("Can only compute Householder matrix for column vectors");
		const alpha = mat.get(0, 0);
		let scale = (mat.rows === 1) ? 0 : Matrix.columnNorm(mat.submatrix(1, mat.rows, 0, 1), 0);
		scale *= scale;
		let v = mat.copy();

		let tau;
		if (scale === 0) {
			tau = 0;			
		} else {
			const t = Math.sqrt(alpha * alpha + scale);
			v.mat[0][0] = (alpha <= 0) ? alpha - t : -scale / (alpha + t);
			tau = 2 * v.get(0, 0) * v.get(0, 0) / (scale + v.get(0, 0) * v.get(0, 0));
			v = Matrix.scale(v, 1 / v.get(0, 0));
		}

		return Matrix.sub(Matrix.identity(mat.rows), Matrix.scale(Matrix.multiply(v, Matrix.transpose(v)), tau));
	}

	static QRDecomposition(mat) {
		/**
		 * Calculates the QR decomposition (https://en.wikipedia.org/wiki/QR_decomposition)
		 * Uses householder reflections
		 */
		let R = mat.copy();
		let Q = Matrix.identity(mat.rows);

		for (let j=0; j<mat.cols; j++) {
			const reflect = Matrix.householderReflection(R.submatrix(j, mat.rows, j, j + 1));
			const H = Matrix.identity(mat.rows);
			for (let a=j; a<mat.rows; a++) {
				for (let b=j; b<mat.rows; b++) {
					H.mat[a][b] = reflect.get(a - j, b - j);
				}
			}
			R = Matrix.multiply(H, R);
			Q = Matrix.multiply(H, Q);
		}

		return {
			Q: Matrix.transpose(Q.submatrix(0, mat.cols, 0, mat.cols)),
			R: R.submatrix(0, mat.cols, 0, mat.cols),
		};
	}

	static eigenvalues(mat) {
		/**
		 * compute eigenvalues of the matrix
		 */
		if (!mat.isSquare) throw new Error("can't compute eigenvalues of non-square matrix");
		let A = mat.copy();
		for (let i=0; i<10; i++) {
			// 10 iterations because I don't want to worry too much about convergence and add so much complexity
			const decomposition = Matrix.QRDecomposition(A);
			A = Matrix.multiply(decomposition.R, decomposition.Q);
		}
		const result = [];
		for (let i=0; i<mat.rows; i++) {
			result.push(A.get(i, i));
		}
		return result;
	}

	static eigenpairs(mat) {
		/** computes the eigenvalues and eigenvectors of the given matrix */
	}

	static exp(mat) {
		/**
		 * computes e to the power of mat
		 * defined by the power series for exp(x)
		 * exp(A) = \sum_{k=0}^{\infty} A^k/k!
		 */
	}

	static log(mat) {
		/**
		 * computes the matrix logarithm
		 * https://en.wikipedia.org/wiki/Logarithm_of_a_matrix
		 */
	}

	static pow(mat) {
		/**
		 * computes the matrix exponential
		 */
	}

	static LUDecomposition(mat) {
		// get the L and U triangular matrices where mat=LU
		// (this implements the Doolittle algorithm)
		if (mat.isSquare) {
			let L, U;
			L = Matrix.emptyMatrix(mat.rows, mat.cols);
			U = Matrix.emptyMatrix(mat.rows, mat.cols);

			for (let i=0; i<mat.rows; i++) {
				for (let k=i; k<mat.rows; k++) {
					// Upper triangular matrix
					let sum = 0;
					for (let j=0; j<i; j++) { sum += L.mat[i][j] * U.mat[j][k]; }
					U.mat[i][k] = mat.mat[i][k] - sum;
				}

				for (let k=i; k<mat.rows; k++) {
					// Lower triangular matrix
					if (i === k) {
						L.mat[i][i] = 1;
					} else {
						let sum = 0;
						for (let j=0; j<i; j++) { sum += L.mat[k][j] * U.mat[j][i]; }
						L.mat[k][i] = (mat.mat[k][i] - sum) / U.mat[i][i];
					}
				}
			}
			return {"L": L, "U": U,};
		}
	throw new Error("The LU-decomposition can only be found for a square matrix.");
	}

	static minors(mat) {
		// return the matrix of minors for mat
		return mat.minors();
	}

	static cofactors(mat) {
		// return the matrix of cofactors for mat
		return mat.cofactors();
	}

	static adjugate(mat) {
		// return the adjugate matrix for mat
		return mat.adjugate();
	}

	static inverse(mat) {
		// return the inverse matrix of mat
		return mat.inverse();
	}

	static transpose(mat) {
		// return the transpose of mat
		return mat.transpose();
	}

	static toLatex(mat) {
		return mat.toLatex();
	}

	static cheapMatrixRound(mat) {
		// replace entries that are within EPSILON of an integer with an integer
		let newMat = mat.copy();
		newMat.applyIndexFunction((v, i, j) => cheapNumberRound(v));
		return newMat;
	}

	copy() {
		// return a copy of this matrix
        const newMatrix = Matrix.emptyMatrix(this.rows, this.cols);
        for (let i=0; i<this.rows; i++) {
            for (let j=0; j<this.cols; j++) {
                newMatrix.mat[i][j] = this.mat[i][j];
            }
        }
		return newMatrix;
	}

	scale(k) {
		// multiply every element in the matrix by k
		for (let i=0; i<this.rows; i++) {
			for (let j=0; j<this.cols; j++) {
				this.mat[i][j] *= k;
			}
		}
	}

	add(mat) {
		if (this.rows === mat.rows && this.cols == mat.cols) {
			for (let i=0; i<this.rows; i++) {
				for (let j=0; j<this.cols; j++) {
					this.mat[i][j] += mat.mat[i][j];
				}
			}
			return;
		}
		throw new Error("Invalid matrix dimensions (both matrices should be the same size).");
	}

	sub(mat) {
		let matCopy = mat.copy();
		matCopy.scale(-1);
		this.add(matCopy);
	}

	applyIndexFunction(f) {
		// replace this.mat[i][j] with f(this.mat[i][j], i, j)
		for (let i=0; i<this.rows; i++) {
			for (let j=0; j<this.cols; j++) {
				this.mat[i][j] = f(this.mat[i][j], i, j);
			}
		}
	}

	minors() {
		// compute the matrix of minors for the matrix
		let res = Matrix.emptyMatrix(this.rows, this.cols);
		for (let i=0; i<this.rows; i++) {
			for (let j=0; j<this.cols; j++) {
				res.mat[i][j] = this.withoutRow(i).withoutColumn(j).determinant();
			}
		}
		return Matrix.cheapMatrixRound(res);
	}

	cofactors() {
		// compute the matrix of cofactors for the matrix
		let minors = this.minors();
		minors.applyIndexFunction((v, i, j) => {
			if ((i+j)%2 === 0) {
				return v;
			} else {
				return -v;
			}
		});
		// could do (v,i,j)=>v*Math.pow(-1,(i+j)%2!==0) but nah
		return minors;
	}

	adjugate() {
		// compute the adjugate/adjoint matrix
		return this.cofactors().transpose();
	}

	inverse() {
		// compute the inverse of the matrix, if possible
		//  (1/det) * adjugate
		let det = this.determinant(); // this will throw if !this.isSquare
		if (det === 0) {
			throw new Error("Matrix is singular.");
		} else {
			return Matrix.cheapMatrixRound(Matrix.scale(this.adjugate(), 1/det));
		}
	}

	transpose() {
		let mat = Matrix.emptyMatrix(this.cols, this.rows);
		for (let i=0; i<this.rows; i++) {
			for (let j=0; j<this.cols; j++) {
				mat.mat[j][i] = this.mat[i][j];
			}
		}
		return mat;
	}

	getColumn(i) {
		// return the ith column
		let col = [];
		for (let j=0; j<this.mat.length; j++) {
			col.push(this.mat[j][i]);
		}
		return col;
	}

	getRow(i) {
		// return the ith row
		return this.mat[i];
	}

	get(i, j) {
		/** Return the element in the ith row and jth column */
		return this.mat[i][j];
	}

	rowSwap(i1, i2) {
		// swap two rows in the matrix
		if (i1 < this.rows && i2 < this.rows) {
			let temp;
			temp = this.mat[i1].slice();
			this.mat[i1] = this.mat[i2].slice();
			this.mat[i2] = temp;
			return;
		}
		throw new Error("Row is invalid.");
	}

	columnSwap(j1, j2) {
		// swap two columns in the matrix
		if (j1 < this.cols && j2 < this.cols) {
			let mat = this.transpose();
			mat.rowSwap(j1, j2);
			this.mat = mat.transpose().mat;
			return;
		}
		throw new Error("Column is invalid.");
	}

	submatrix(i1, i2, j1, j2) {
		/**
		 * return the submatrix [ a_ij ] with i1 <= i < i2, j1 <= j < j2
		 * */
		if (!(0 <= i1 && i1 <= this.rows && 0 <= j1 && j1 <= this.cols)) throw new Error("submatrix index out of range");
		if (i1 > i2 || j1 > j2) throw new Error("invalid submatrix indices");

		const result = [];
		for (let i=i1; i<i2; i++) {
			const row = [];
			for (let j=j1; j<j2; j++) {
				row.push(this.mat[i][j]);
			}
			result.push(row);
		}
		return new Matrix(result);
	}

	withoutRow(i) {
		// return a copy of the matrix without it's ith row
		if (i >= this.rows) throw new Error("Row is invalid.");
		let mat = this.copy().mat;
		mat.splice(i, 1);
		return new Matrix(mat);
	}

	withoutColumn(j) {
		if (j >= this.cols) throw new Error("Column is invalid.");
		let mat = this.transpose();
		return mat.withoutRow(j).transpose();
	}

	scaleRow(i, k) {
		if (i < this.rows) {
			for (let j=0; j<this.cols; j++) {
				this.mat[i][j] *= k;
			}
			return;
		}
		throw new Error("Invalid row.");
	}

	subtractScaledRow(row, rowToBeSubtracted, scale) {
		// subtract scale * rowToBeSubtracted from row
		rowToBeSubtracted = this.mat[rowToBeSubtracted];
		for (let j=0; j<this.cols; j++) {
			this.mat[row][j] -= scale * rowToBeSubtracted[j];
		}
	}

	rowHasLeadingOne(i) {
		// return if the column has a leading one
		if (i < this.rows) {
			for (let j=0; j<this.cols; j++) {
				if (this.mat[i][j] === 1) {
					return true;
				} else if (this.mat[i][j] !== 0) {
					return false;
				}
			}
		}
		throw new Error("Invalid row.");
	}

	columnOfFirstNonzero(i) {
		// return the column of the first nonzero entry in the ith row (if none: return -1)
		if (i < this.rows) {
			for (let j=0; j<this.cols; j++) {
				if (this.mat[i][j] !== 0) return j;
			}
			return -1;
		}
		throw new Error("Invalid row.");
	}

	ref() {
		// convert the matrix to row-echelon form

	}

	rref() {
		// convert the matrix to reduced row-echelon form
		//this.ref();
		let lead = 0, i;
		for (let r=0; r<this.rows; r++) {
			if (this.cols <= lead) break;
			i = r;
			while (this.mat[i][lead] === 0) {
				i++;
				if (this.rows === i) {
					i = r;
					lead++;
					if (this.cols === lead) break;
				}
			}
			this.rowSwap(i, r);
			if (this.mat[r][lead] !== 0) this.scaleRow(r, 1/this.mat[r][lead]);
			for (i=0; i<this.rows; i++) {
				if (i !== r) this.subtractScaledRow(i, r, this.mat[i][lead]);
			}
			lead++;
		}
	}

	determinant() {
		// compute the determinant of the matrix
		return Matrix.det(this);
	}

	trace() {
		// compute the trace of the matrix
		if (this.isSquare) {
			let total = 0;
			for (let i=0; i<this.rows; i++) {
				total += this.mat[i][i];
			}
			return total;
		}
		throw new Error("Trace is only defined for square matrices.");
	}

    rowSum(i) {
        // get the sum of the ith row
        let total = 0;
        const row = this.getRow(i);
        return row.reduce((x, y) => x + y, 0);
    }

    columnSum(j) {
        // get the sum of the jth column
        let total = 0;
        const col = this.getColumn(j)
        return col.reduce((x, y) => x + y, 0);
    }

    rowAbsoluteSum(i) {
        // get the sum of the absolute values of the ith row
        let total = 0;
        const row = this.getRow(i);
        return row.reduce((x, y) => Math.abs(x) + Math.abs(y), 0);
    }

    columnAbsoluteSum(j) {
        // get the sum of the absolute values of the jth column
        let total = 0;
        const col = this.getColumn(j)
        return col.reduce((x, y) => Math.abs(x) + Math.abs(y), 0);
    }

    rowSums() {
        // get the sum of every row
        const totals = [];
        for (let i=0; i<this.rows; i++) {
            totals.push(this.rowSum(i));
        }
        return totals;
    }

    columnSums() {
        // get the sum of every column
        const totals = [];
        for (let j=0; j<this.cols; j++) {
            totals.push(this.columnSum(j));
        }
        return totals;
    }

    rowAbsoluteSums() {
        // get the absolute sum of every row
        const totals = [];
        for (let i=0; i<this.rows; i++) {
            totals.push(this.rowAbsoluteSum(i));
        }
        return totals;
    }

    columnAbsoluteSums() {
        // get the absolute sum of every column
        const totals = [];
        for (let j=0; j<this.cols; j++) {
            totals.push(this.columnAbsoluteSum(j));
        }
        return totals;
    }

    norm() {
        /*
        compute the 1-norm of the matrix
        |A| = sup_{x of norm 1} |Ax| = max_{x of norm 1} |Ax| (since R^n is proper)
        */
       return Math.max(...this.columnAbsoluteSums());
    }

	toString() {
		// return the string representation of the matrix
		let str = "", seg;
		for (let i=0; i<this.rows; i++) {
			str += "|  ";
			for (let j=0; j<this.cols; j++) {
				seg = "" + Math.floor(this.mat[i][j] * 100) / 100;
				seg += " ".repeat(5-seg.length);
				str += seg;
				if (j < this.cols - 1) str += "   ";
			}
			str += "  |\n"
		}
		str = " _" + " ".repeat(3 * (this.cols - 1) + 5 * this.cols+2) + "_\n" + str;
		str += " -" + " ".repeat(3 * (this.cols - 1) + 5 * this.cols+2) + "- ";
		return str;
	}

	toLatex(openDelim="$$", closeDelim="$$") {
		// return the pmatrix form of the matrix
		let str = openDelim + "\\begin{pmatrix}";
		for (let i=0; i<this.rows; i++) {
			for (let j=0; j<this.cols; j++) {
				if (j+1 < this.cols) {
					str += this.mat[i][j] + "&"
				} else {
					str += "" + this.mat[i][j];
				}
			}
			if (i !== this.rows-1) str += "\\\\";
		}
		str += "\\end{pmatrix}" + closeDelim;
		return str;
	}

	show() {
		// this relies on the fact that the default inspect window in most browsers uses a monospaced font
		let str = this.toString();
		console.log(str);
	}

}

function matrix(arr) {
    // utility function for easy matrix creation
	return new Matrix(arr);
}


module.exports = {
	Matrix, matrix,
};

},{}],7:[function(require,module,exports){


const { matrix } = require("./matrix.js");
const { complex } = require("./complex.js");


function stereographic(z) {
    /**
     * computes the stereographic projection of the 3d point z onto the xy plane
     */
    return complex(
        z[0] / (1 - z[2]),
        z[1] / (1 - z[2]),
    );
}

function inverseStereoProject(z) {
    /**
     * computes the inverse stereographic projection of a 2d point onto the unit spheree
     */
    const normSq = z.normSq();
        return matrix([
            (2 * z.re) / (1 + normSq),
            (2 * z.im) / (1 + normSq),
            (normSq - 1) / (1 + normSq),
        ]).transpose();
}

function perspectiveProject(Z, orbit, fov) {
    /**
     * computes the perspective projection of the 3d point Z onto the xz plane
     * given the orbit radius and fov values
     */
    const values = Z.getColumn(0);
    const x = values[0], y = values[1], z = values[2];
    const denom = (orbit + fov * y);
    return complex(
        x / denom, z / denom,
    );
}


module.exports = {
    stereographic, inverseStereoProject,
    perspectiveProject,
};
},{"./complex.js":2,"./matrix.js":6}],8:[function(require,module,exports){



class rVector {

    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    dot(v) {
        return this.x * v.x + this.y * v.y;
    }

    magSq() {
        return this.x * this.x + this.y * this.y;
    }

    mag() {
        return Math.sqrt(this.magSq());
    }

    scale(k) {
        return new rVector(this.x * k, this.y * k);
    }

    proj(v) {
        return v.scale(this.dot(v) / v.magSq());
    }

    unit() {
        return this.scale(1 / this.mag());
    }

    angleBetween(v) {
        return Math.acos(this.dot(v).scale(this.mag() * v.mag()));
    }

    perp() {
        return new rVector(-this.y, this.x);
    }

    reflect(pointOnMirror, mirrorNormal) {
        return pointOnMirror.add(this.sub(pointOnMirror).proj(mirrorNormal)).scale(2).sub(this);
    }

    add(v) {
        return new rVector(this.x + v.x, this.y + v.y);
    }

    sub(v) {
        return new rVector(this.x - v.x, this.y - v.y);
    }

}


function rvec(x, y) {
    return new rVector(x, y);
}


module.exports = {
    rVector, rvec,
};
},{}],9:[function(require,module,exports){



class ErrorTracker {

    constructor(errorDivID, successMsg, callback=null, successCallback=null) {
        this.hasError = false;
        this.message = null;
        this.successMsg = successMsg;
        this.target = errorDivID;
        this.callback = callback;
        this.successCallback = successCallback;
    }

    setTarget(errorDivID) {
        this.target = errorDivID;
    }

    setCallback(callback) {
        this.callback = callback;
    }

    setSuccessCallback(successCallback) {
        this.successCallback = successCallback;
    }

    defaultCallback(target) {
        if (target) {
            target.innerText = this.message;
            target.style.color = "red";
            target.style.display = "block";
        }
    }

    defaultSuccessCallback(target) {
        if (target) {
            target.innerText = this.message;
            target.style.color = "green";
            target.style.display = "block";
        }
    }

    error(message) {
        this.hasError = true;
        this.message = message;

        const target = document.querySelector(`#${this.target}`);

        if (this.callback === null) {
            this.defaultCallback(target);
        } else {
            this.callback(this.message, target);
        }
    }

    clear() {
        this.hasError = false;
        this.message = this.successMsg;

        const target = document.querySelector(`#${this.target}`);
        
        if (this.successCallback === null) {
            this.defaultCallback(target);
        } else {
            this.successCallback(this.message, target);
        }
    }

}


const tracker = new ErrorTracker(null, "Parsing successful!");

module.exports = {
    tracker,
};
},{}],10:[function(require,module,exports){
function clearFractions(text, depth=0) {
    /*
    convert \frac{num}{denom} to ((num)/(denom))
    */
    if (depth > 9) {
        console.error("fractions nested too deeply");
        return "";
    }
    let newString = "", letter, buffer, bCount, j, fracOpen=false;
    for (let i=0; i<text.length; i++) {
        letter = text[i];
        if (i >= 4 && text.slice(i-4, i) === "frac") fracOpen = true;
        if (letter === "{" && fracOpen) {
            buffer = "";
            bCount = 1;
            for (j=i+1; j<text.length; j++) {
                if (text[j] === "{") {
                    bCount++;
                } else if (text[j] === "}") {
                    bCount--;
                if (bCount === 0) break;
                }
                buffer += text[j];
            }
            if (i >= 4 && text.slice(i-4, i) === "frac") {
                newString += "((" + clearFractions(buffer, depth+1) +")/"
            } else {
                newString += "(" + clearFractions(buffer, depth+1) + "))";
                fracOpen = false;
            }
            i = j;
        } else {
            newString += letter;
        }
    }
    return newString.replace(/\\frac/g, "");
}

/*
sum and prod will be written as functions e.g.
sum(k, [1,...,infinity], 1/k^2)
prod(p, [2, 3, 5, 7, 11, 13], 1/(1-1/p))
*/
function cleanLatex(text) {
    text = text.replace(/\s/g, "");                                                                 // remove whitespace
    text = text.replace(/\\cdot/g, "*").replace(/\\operatorname{(.*?)}/g, (tot, group1)=>group1);   // remove operator wrappers
    text = clearFractions(text);                                                                    // convert latex fractions to a/b
    text = text.replace(/(\\left)|(\\right)/g, "");                                                 // remove \left and \right
    text = text.replace(/\\/g, "");                                                                 // remove backslashes
    text = text.replace(/{/g, "(").replace(/}/g, ")");                                              // replace curly brackets with parentheses
    return text;
}


/*
test strings

4\cos\left(a+bi\right)\cdot\sin\left(x\right)+\frac{1}{2iz}
\frac{1}{1+\frac{1}{x}+\frac{1}{y}}
*/



module.exports = {
    cleanLatex,
};
},{}],11:[function(require,module,exports){

const { Parser } = require("./parser.js");
const {
    InfixParselet, PrefixParselet,
    AssignParselet, BinaryOperatorParselet,
    CallParselet, GroupParselet,
    NameParselet, NumberParselet,
    PrefixOperatorParselet,
} = require("./parselets.js");
const { TokenType } = require("./tokentype.js");
const { Precedence } = require("./precedence.js");



class ExpressionParser extends Parser {

    constructor(tokens) {
        super(tokens);

        this.registerPrefix(TokenType.NUMBER, new NumberParselet());
        this.registerPrefix(TokenType.NAME, new NameParselet());

        this.registerInfix(TokenType.ASSIGN, new AssignParselet());
        this.registerPrefix(TokenType.LEFT_PAREN, new GroupParselet());
        this.registerInfix(TokenType.LEFT_PAREN, new CallParselet());

        this.prefix(TokenType.PLUS, Precedence.PREFIX);
        this.prefix(TokenType.MINUS, Precedence.PREFIX);

        this.infixLeft(TokenType.COLON, Precedence.CALL);
        this.infixLeft(TokenType.PLUS, Precedence.SUM);
        this.infixLeft(TokenType.MINUS, Precedence.SUM);
        this.infixLeft(TokenType.ASTERISK, Precedence.PRODUCT);
        this.infixLeft(TokenType.SLASH, Precedence.PRODUCT);
        this.infixRight(TokenType.CARET, Precedence.EXPONENT);

    }

    prefix(tokenType, precedence) {
        this.registerPrefix(tokenType, new PrefixOperatorParselet(precedence));
    }

    infixLeft(tokenType, precedence) {
        this.registerInfix(tokenType, new BinaryOperatorParselet(precedence, false));
    }

    infixRight(tokenType, precedence) {
        this.registerInfix(tokenType, new BinaryOperatorParselet(precedence, true));
    }

}


module.exports = {
    ExpressionParser,
};
},{"./parselets.js":14,"./parser.js":15,"./precedence.js":16,"./tokentype.js":18}],12:[function(require,module,exports){

const { TokenType } = require("./tokentype.js");



class Expression {

}


class AssignExpression extends Expression {

    constructor(left, right) {
        super();
        this.mLeft = left;
        this.mRight = right;
    }

    toString() {
        return `(${this.mLeft.toString()} = ${this.mRight.toString()})`
    }

}


class CallExpression extends Expression {

    constructor(func, args) {
        super();
        this.mFunction = func;
        this.mArgs = args;
    }

    toString() {
        const args = this.mArgs.join();
        return `${this.mFunction.toString()}(${args})`;
    }

}


class NameExpression extends Expression {

    constructor(name) {
        super();
        this.mName = name;
    }

    toString() {
        return this.mName;
    }

}


class NumberExpression extends Expression {

    constructor(number) {
        super();
        this.mNumber = number;
    }

    toString() {
        return this.mNumber.toString();
    }

}


class OperatorExpression extends Expression {
    
    constructor(left, operator, right) {
        super();
        this.mLeft = left;
        this.mOperator = operator;
        this.mRight = right;
    }

    toString() {
        return `(${this.mLeft.toString()} ${TokenType.toStr(this.mOperator)} ${this.mRight.toString()})`
    }

}


class PrefixExpression extends Expression {

    constructor(operator, right) {
        super();
        this.mOperator = operator;
        this.mRight = right;
    }

    toString() {
        return `(${TokenType.toStr(this.mOperator)}${this.mRight})`;
    }

}


module.exports = {
    Expression,
    AssignExpression,
    CallExpression,
    NameExpression,
    NumberExpression,
    OperatorExpression,
    PrefixExpression,
};
},{"./tokentype.js":18}],13:[function(require,module,exports){

const { TokenType } = require("./tokentype.js");
const { Token } = require("./token.js");
const { Trie } = require("../trie.js");
const { tracker } = require("../errors.js");


const num = "0123456789";
const alnum = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";



class Lexer {

    constructor(mText, allowUnboundIdentifiers=false) {
        this.mPunctuators = [];
        this.mText = mText;
        this.index = 0;
        this.allowUnboundIdentifiers = allowUnboundIdentifiers;
        this.setScope({});
    }

    setScope(scope) {
        this.scope = scope;
        this.builtinLookup = new Trie(Object.keys(scope.builtin === undefined ? {} : scope.builtin));
        this.userGlobalLookup = new Trie(Object.keys(scope.userGlobal === undefined ? {} : scope.userGlobal));
    }

    setAllowUnboundIdentifiers(allowUnboundIdentifiers) {
        this.allowUnboundIdentifiers = allowUnboundIdentifiers;
    }

    setText(mText) {
        this.mPunctuators = [];
        this.mText = mText;
        this.index = 0;
    }

    getTokenName(tokenizingAssignment=false, assignmentEncountered=false) {
        // consume characters until non-identifier character is encountered, add them to buffer
        let buffer = "";
        while (this.index < this.mText.length && alnum.includes(this.mText[this.index])) {
            buffer += this.mText[this.index];
            this.index++;
        }

        // split buffer into identifiers, add implicit multiplication where necessary
        const identifiers = [];
        while (buffer.length > 0) {
            let matchFound = false;
            let possibleIdentifier, i;
            for (i=buffer.length-1; i>=0; i--) {
                possibleIdentifier = buffer.slice(0, i+1);
                if (this.builtinLookup.containsKey(possibleIdentifier) || this.userGlobalLookup.containsKey(possibleIdentifier)) {
                    matchFound = true;
                    break;
                }
            }
            if (matchFound) {
                identifiers.push(possibleIdentifier);
                buffer = buffer.slice(i+1, buffer.length);
            } else {
                if (alnum.includes(possibleIdentifier[0]) && !num.includes(possibleIdentifier[0])) {
                    if (this.allowUnboundIdentifiers || (tokenizingAssignment)) {
                        // no match found in the remaining part of the buffer, so greedily make the 
                        // unidentified characters a single identifier
                        identifiers.push(buffer);
                        buffer = "";
                    } else {
                        tracker.error(`Undefined identifier ${buffer}`);
                        this.kill();
                        return;
                    }
                } else {
                    tracker.error("unexpected number following identifier");
                    this.kill();
                    return;
                }
            }
        }
        for (let identifier of identifiers) {
            this.mPunctuators.push(new Token(TokenType.NAME, identifier));
            this.mPunctuators.push(new Token(TokenType.ASTERISK, "*"));
        }
        this.mPunctuators.pop();
    }

    getTokenNumber() {
        let number = "";
        let decimalEncounted = false;
        while (this.index < this.mText.length && (num.includes(this.mText[this.index]) || this.mText[this.index] === ".")) {
            if (this.mText[this.index] === ".") {
                if (decimalEncounted) {
                    tracker.error("two decimal points encountered in number");
                    this.kill();
                    return;
                }
                decimalEncounted = true;
            }
            number += this.mText[this.index];
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.NUMBER, number));
    }

    kill() {
        /** make tokenization stop */
        this.index = this.mText.length;
    }

    tokenize() {
        this.mPunctuators = [];
        this.index = 0;
        const tokenizingAssignment = this.mText.includes("=");
        let assignmentEncountered = false;

        while (this.index < this.mText.length) {
            const char = this.mText[this.index];
            
            if (char === "(") {
                this.mPunctuators.push(new Token(TokenType.LEFT_PAREN, char));
            } else if (char === ")") {
                this.mPunctuators.push(new Token(TokenType.RIGHT_PAREN, char));
            } else if (char === ",") {
                this.mPunctuators.push(new Token(TokenType.COMMA, char));
            } else if (char === "=") {
                this.mPunctuators.push(new Token(TokenType.ASSIGN, char));
                assignmentEncountered = true;
            } else if (char === "+") {
                this.mPunctuators.push(new Token(TokenType.PLUS, char));
            } else if (char === "-") {
                this.mPunctuators.push(new Token(TokenType.MINUS, char));
            } else if (char === "*") {
                this.mPunctuators.push(new Token(TokenType.ASTERISK, char));
            } else if (char === "/") {
                this.mPunctuators.push(new Token(TokenType.SLASH, char));
            } else if (char === "^") {
                this.mPunctuators.push(new Token(TokenType.CARET, char));
            } else if (char === ":") {
                this.mPunctuators.push(new Token(TokenType.COLON, char));
            } else if (num.includes(char) || char === ".") {
                this.getTokenNumber();
                continue;
            } else if (alnum.includes(char)) {
                this.getTokenName(tokenizingAssignment, assignmentEncountered);
                continue;
            } else {
                tracker.error(`bruh what is this character even doing in your input ${char}`);
                this.kill();
                return;
            }
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.EOF, null));

        /**
         * Implicit multiplication loop
         * cases handled here (let x, y = identifiers and n = number ):
         * xn -> syntax error
         * nx -> n * x
         * n( -> n * (
         * )n -> ) * n
         * )x -> ) * x
         * )( -> ) * (
         * x( -> x(         (assuming scope[x].isFunction === false)
         * case not handled here:
         * xy -> x * y. this case is handled in getTokenName()
         */
        const tokens = (this.mPunctuators.length > 0) ? [this.mPunctuators[0]] : [];
        for (let i=0; i<this.mPunctuators.length-1; i++) {
            const token = this.mPunctuators[i];
            const nextToken = this.mPunctuators[i+1];
            let needsMultiplication = false;

            if (token.mtype === TokenType.NAME && nextToken.mtype === TokenType.NUMBER) {
                // this case will theoretically never happen (it'll raise an error in getTokenName()), but it's here for completeness
                tracker.error("unexpected number following identifier");
                this.kill();
                return;
            } else if ( // obviously, you can reduce the number of comparisons in this condition, but it's more clear this way
                token.mtype === TokenType.NUMBER && nextToken.mtype === TokenType.NAME ||           // nx
                token.mtype === TokenType.NUMBER && nextToken.mtype === TokenType.LEFT_PAREN ||     // n(
                token.mtype === TokenType.RIGHT_PAREN && nextToken.mtype === TokenType.NUMBER ||    // )n
                token.mtype === TokenType.RIGHT_PAREN && nextToken.mtype === TokenType.NAME ||      // )x
                token.mtype === TokenType.RIGHT_PAREN && nextToken.mtype === TokenType.LEFT_PAREN   // )(
            ) {
                needsMultiplication = true;
            } else if (token.mtype === TokenType.NAME && nextToken.mtype === TokenType.LEFT_PAREN) {
                if (this.builtinLookup.containsKey(token.text)) {
                    needsMultiplication = !this.scope.builtin[token.text].isFunction;
                } else if (this.userGlobalLookup.containsKey(token.text)) {
                    needsMultiplication = !this.scope.userGlobal[token.text].isFunction;
                }
            }

            if (needsMultiplication) tokens.push(new Token(TokenType.ASTERISK, "*"));
            tokens.push(nextToken);
        }

        this.mPunctuators = tokens.slice();
        this._validateTokens();
    }

    _checkIsFunction(token) {
        return (
            token.mtype === TokenType.NAME && 
            (!!this.scope.builtin[token.text]?.isFunction || !!this.scope.userGlobal[token.text]?.isFunction)
        );
    }

    _validateTokens() {
        let equalCount = 0;
        for (let i=0; i<this.mPunctuators.length-1; i++) {
            /**
             * in the case that the input string is empty, this.mPunctuators contains only the EOL token
             * in all other cases, this loop will run at least once
             */
            const token = this.mPunctuators[i];
            const nextToken = this.mPunctuators[i+1];

            if (this._checkIsFunction(token) && nextToken.mtype !== TokenType.LEFT_PAREN) {
                tracker.error(`Function call to ${token.text} requires parentheses`);
                return;
            }

            if (token.mtype === TokenType.ASSIGN) equalCount++;
        }

        if (equalCount > 1) {
            tracker.error("More than one equals sign present in expression");
        }
    }

    getTokens() {
        return this.mPunctuators;
    }

}


module.exports = {
    Lexer,
};
},{"../errors.js":9,"../trie.js":19,"./token.js":17,"./tokentype.js":18}],14:[function(require,module,exports){


const { Precedence } = require("./precedence.js");
const { TokenType } = require("./tokentype.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("./expressions.js");
const { tracker } = require("../errors.js");


class InfixParselet {

    parse(parser, left, token) {

    }

    getPrecedence() {

    }

}


class PrefixParselet {

    parse(parser, token) {

    }

    getPrecedence() {

    }

}


class AssignParselet extends InfixParselet {

    parse(parser, left, token) {
        const right = parser.parseExpression(Precedence.ASSIGNMENT - 1);
        if (!(left instanceof NameExpression || left instanceof CallExpression)) {
            tracker.error("lhs of assignment should be identifier or function name with arguments");
            return;
        }

        return new AssignExpression(left, right);
    }

    getPrecedence() {
        return Precedence.ASSIGNMENT;
    }

}


class BinaryOperatorParselet extends InfixParselet {

    constructor(precedence, mIsRight) {
        super();
        this.mPrecedence = precedence;
        this.mIsRight = mIsRight;
    }

    parse(parser, left, token) {
        const right = parser.parseExpression(this.mPrecedence - (this.mIsRight ? 1 : 0));
        return new OperatorExpression(left, token.mtype, right);
    }

    getPrecedence() {
        return this.mPrecedence;
    }

}


class CallParselet extends InfixParselet {

    constructor() {
        super();
    }

    parse(parser, left, token) {
        const args = [];
        
        if (!parser.peek(TokenType.RIGHT_PAREN)) {
            while (!tracker.hasError) {
                args.push(parser.parseExpression());
                if (parser.peek(TokenType.RIGHT_PAREN)) {
                    parser.consume(TokenType.RIGHT_PAREN);
                    break;
                }
                parser.consume(TokenType.COMMA);
            }
        }

        return new CallExpression(left, args);
    }

    getPrecedence() {
        return Precedence.CALL;
    }

}


class GroupParselet extends PrefixParselet {

    parse(parser, token) {
        const expr = parser.parseExpression();
        parser.consume(TokenType.RIGHT_PAREN);
        return expr;
    }

}


class NameParselet extends PrefixParselet {

    constructor() {
        super();
    }

    parse(parser, token) {
        return new NameExpression(token.text);
    }

}


class NumberParselet extends PrefixParselet {

    constructor() {
        super();
    }

    parse(parser, token) {
        const value = parseFloat(token.text);
        if (isNaN(value)) {
            tracker.error(`Invalid number literal ${value}`);
        } else {
            return new NumberExpression(parseFloat(token.text));
        }
    }

}


class PrefixOperatorParselet extends PrefixParselet {

    constructor(precedence) {
        super();
        this.mPrecedence = precedence;
    }

    parse(parser, token) {
        const right = parser.parseExpression(this.mPrecedence);
        return new PrefixExpression(token.mtype, right);
    }

}



module.exports = {
    InfixParselet, PrefixParselet,
    AssignParselet, BinaryOperatorParselet,
    CallParselet, GroupParselet,
    NameParselet, NumberParselet,
    PrefixOperatorParselet,
};
},{"../errors.js":9,"./expressions.js":12,"./precedence.js":16,"./tokentype.js":18}],15:[function(require,module,exports){



const { tracker } = require("../errors.js");
const { TokenType } = require("./tokentype.js");
const { Precedence } = require("./precedence.js");
const { Token } = require("./token.js");


class Parser {

    constructor(tokens) {
        this.mTokens = tokens;
        this.mPrefixParselets = {};
        this.mInfixParselets = {};
        this.index = 0;
    }

    registerPrefix(tokenType, parselet) {
        this.mPrefixParselets[tokenType] = parselet;
    }

    registerInfix(tokenType, parselet) {
        this.mInfixParselets[tokenType] = parselet;
    }

    parseExpression(precedence=Precedence.LOWEST) {
        let token = this.consume();
        if (tracker.hasError) return;
        if (!Object.keys(this.mPrefixParselets).includes(token.mtype.toString())) {
            if (token.mtype === TokenType.EOF) {
                tracker.error(`Unexpected EOF while parsing (1)`);
            } else {
                tracker.error(`Invalid syntax`);
            }
            return;
        }

        let left = this.mPrefixParselets[token.mtype].parse(this, token);
        if (tracker.hasError) return;

        while (precedence < this.getPrecedence()) {
            token = this.consume();
            if (Object.keys(this.mInfixParselets).includes(token.mtype.toString())) {
                left = this.mInfixParselets[token.mtype].parse(this, left, token);
                if (tracker.hasError) return;
            }
        }

        return left;
    }

    consume(expected=null) {
        const token = this.lookAhead(0);
        if (expected !== null && token.mtype === TokenType.EOF) {
            tracker.error(`Unexpected EOF while parsing`);
            return;
        } else if (expected !== null && token.mtype !== expected) {
            tracker.error(`Unexpected EOF while parsing`);
            return;
        }
        this.index++;
        return token;
    }

    lookAhead(distance) {
        if (this.index + distance < this.mTokens.length) {
            return this.mTokens[this.index + distance];
        } else {
            return new Token(TokenType.EOF, null);
        }
    }

    peek(expected) {
        const token = this.lookAhead(0);
        return (token.mtype === expected);
    }

    getPrecedence() {
        const currentToken = this.lookAhead(0);
        if (Object.keys(this.mInfixParselets).includes(currentToken.mtype.toString())) {
            return this.mInfixParselets[currentToken.mtype].getPrecedence();
        }
        return Precedence.LOWEST;
    }

}


module.exports = {
    Parser,
};
},{"../errors.js":9,"./precedence.js":16,"./token.js":17,"./tokentype.js":18}],16:[function(require,module,exports){



class Precedence {

    static LOWEST = 0;
    static ASSIGNMENT = 1;
    static SUM = 2;
    static PRODUCT = 3;
    static EXPONENT = 4;
    static PREFIX = 5;
    static CALL = 6;

}

module.exports = {
    Precedence,
};
},{}],17:[function(require,module,exports){


class Token {

    constructor(tokenType, text) {
        this.mtype = tokenType;
        this.text = text;
    }

    toString() {
        return `<${this.mtype}, ${this.text}>`;
    }

}


module.exports = {
    Token,
};
},{}],18:[function(require,module,exports){



class TokenType {

    static LEFT_PAREN = 1;
    static RIGHT_PAREN = 2;
    static COMMA = 3;
    static ASSIGN = 4;
    static PLUS = 5;
    static MINUS = 6;
    static ASTERISK = 7;
    static SLASH = 8;
    static CARET = 9;
    static COLON = 10;
    static NAME = 11;
    static NUMBER = 12;
    static EOF = 13;

    constructor(value) {
        this.value = value;
    }

    static toStr(tokenType) {
        const values = [
            "(", ")", ",", "=", "+", "-", "*", "/", "^", ":", "", "", "",
        ];
        return values[tokenType - 1];
    }

}

module.exports = {
    TokenType,
};
},{}],19:[function(require,module,exports){




class Node {

    static nodeID = 0;

    constructor(value) {
        this.value = value;
        this.id = Node.nodeID;
        Node.nodeID++;
        this.next = null;   // horizontal child
        this.child = null;  // vertical child
    }

    setNext(next) {
        this.next = next;
    }

    setChild(child) {
        this.child = child;
    }

    hasNext() {
        return this.next !== null;
    }

    toString() {
        return `Node(${this.value})`;
    }

}


class LinkedList {

    constructor(valuesArray=null) {
        this.head = null;
        this.size = 0;

        if (valuesArray !== null) {
            for (let value of valuesArray) {
                this.push(value);
            }
        }
    }

    push(value) {
        const newNode = new Node(value);
        
        if (this.head === null) {
            this.head = newNode;
        } else {
            let node = this.head;
            while (node.hasNext()) {
                node = node.next;
            }
            node.setNext(newNode);
        }
        this.size++;
        return newNode;
    }

    remove(value) {
        /*
        Remove value from the linked list.
        Return boolean indicating if successfully removed
        */
        if (this.head === null) {
            return false;
        } else {
            if (this.head.value === value) {
                this.head = this.head.next;
                this.size--;
                return true;
            } else {
                let node = this.head;
                while (node.hasNext()) {
                    if (node.next.value === value) {
                        node.next = node.next.next;
                        this.size--;
                        return true;
                    }
                    node = node.next;
                }
                return false;
            }
        }
    }

    search(value) {
        let node = this.head;
        while (node !== null) {
            if (node.value === value) {
                return node;
            }
            node = node.next;
        }
        return null;
    }

    toString() {        
        if (this.head === null) {
            return "LinkedList([])";
        } else {
            let result = "[";
            let node = this.head;
            while (node !== null) {
                result += node.toString() + ", ";
                node = node.next;
            }
            result = result.slice(0, -2) + "]";
            return result;
        }
    }

}


class Trie {

    static TERMINATOR = "<end>";

    constructor(keysArray=null) {
        this.root = null;
        this.maxDepth = 0;
        this.size = 0;

        if (keysArray !== null) {
            for (let key of keysArray) {
                this.insertKey(key);
            }
        }
    }

    insertKey(key) {
        if (this.root === null) {
            this.root = new LinkedList();
        }

        this.maxDepth = Math.max(this.maxDepth, key.length);

        let list = this.root;
        for (let i=0; i<key.length; i++) {
            const character = key[i];
            const node = list.search(character);
            if (node === null) {
                const newNode = list.push(character);
                newNode.setChild(new LinkedList());
                list = newNode.child;
            } else {
                list = node.child;
            }

            if (i === key.length - 1 && list.search(Trie.TERMINATOR) === null) {
                list.push(Trie.TERMINATOR);
                this.size++;
            }
        }
    }

    containsKey(key) {
        if (this.root === null) return false;

        let list = this.root;
        for (let character of key) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        return list.search(Trie.TERMINATOR) !== null;
    }

    containsPrefix(prefix) {
        if (this.root === null) return false;

        let list = this.root;
        for (let character of prefix) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        return true;
    }

    remove(key) {
        /*
        Remove key from trie, return boolean indicating if the removal was successful
        Note: does NOT remove prefixes.
        For example:
        const trie = new Trie(["dog"]);
        trie.remove("dog"); // true
        trie.containsKey("dog"); // false
        trie.containsPrefix("dog"); // true
        */
        if (this.root === null) return false;
        let list = this.root;
        for (let character of key) {
            const node = list.search(character);
            if (node === null) return false;
            list = node.child;
        }
        const removalResult = list.remove(Trie.TERMINATOR);
        if (removalResult) this.size--;
        return removalResult;
    }

    prefixSearch(prefix, maxResults=100) {
        /**
         * Search the trie for keys with the given prefix.
         * Returns no more than maxResults keys
         */
        if (this.root === null || prefix === "") return [];

        let list = this.root, node;
        for (let character of prefix) {
            node = list.search(character);
            if (node === null) return []; // prefix not in trie
            list = node.child;
        }

        const results = [];
        if (list.search(Trie.TERMINATOR) !== null) results.push(prefix); 

        // depth-first search for matching keys
        const queue = [{node: list.head, suffix: ""}];
        let current;
        while (results.length < maxResults && queue.length > 0) {
            current = queue.pop();
            let currentNode = current.node, currentSuffix = current.suffix;
            while (currentNode !== null) {
                if (currentNode.value === Trie.TERMINATOR) {
                    currentNode = currentNode.next;
                    continue;
                }
                const searchResult = currentNode.child.search(Trie.TERMINATOR);
                if (searchResult !== null) results.push(prefix + currentSuffix + currentNode.value);
                queue.unshift({
                    node: currentNode.child.head,
                    suffix: currentSuffix + currentNode.value,
                });
                currentNode = currentNode.next;
                if (results.length === maxResults) break;
            }
        }
        return results;
    }

}

module.exports = {
    Node, LinkedList, Trie,
};
},{}],20:[function(require,module,exports){

const { complex, Complex } = require("./math/complex.js");


const scope = {
    builtin: {
        "z": { isFunction: false, },

        "norm": { isFunction: true, },
        "normSq": { isFunction: true, },
        "arg": { isFunction: true, },
        "inv": { isFunction: true, },
        "exp": { isFunction: true, },
        "ln": { isFunction: true, },
        "sqrt": { isFunction: true, },
        "sin": { isFunction: true, },
        "cos": { isFunction: true, },
        "tan": { isFunction: true, },
        "sinh": { isFunction: true, },
        "cosh": { isFunction: true, },
        "tanh": { isFunction: true, },
        "re": { isFunction: true, },
        "im": { isFunction: true, },
        "Gamma": { isFunction: true, },
        "beta": { isFunction: true, },
        "min": { isFunction: true, },
        "max": { isFunction: true, },
        "lerp": { isFunction: true, },

        "i": { isFunction: false, },
        "pi": { isFunction: false, },
        "tau": { isFunction: false, },
        "e": { isFunction: false, },
        "realBounds": { isFunction: false, },
        "imagBounds": { isFunction: false, },

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
    "sinh": Complex.sinh,
    "cosh": Complex.cosh,
    "tanh": Complex.tanh,
    "re": (z) => complex(z.re, 0),
    "im": (z) => complex(z.im, 0),
    "Gamma": Complex.gamma,
    "beta": Complex.beta,
    "min": Complex.min,
    "max": Complex.max,
    "lerp": (z1, z2, t) => Complex.mult(complex(1, 0).sub(t), z1).add(z2.mult(t)),
  
    "i": complex(0, 1),
    "pi": complex(Math.PI, 0),
    "tau": complex(2 * Math.PI, 0),
    "e": complex(Math.E, 0),
    "realBounds": complex(0, 0),
    "imagBounds": complex(0, 0),
};

const valueScope = {};

module.exports = {
    scope, defaultValueScope, valueScope,
};
},{"./math/complex.js":2}],21:[function(require,module,exports){
/**
 * imports
 */


const { cleanLatex } = require("./parsing/latex_convert.js");
const { tracker } = require("./parsing/errors.js");
const { TokenType } = require("./parsing/pratt/tokentype.js");
const {
    Expression, AssignExpression,
    CallExpression, NameExpression,
    NumberExpression, OperatorExpression,
    PrefixExpression
} = require("./parsing/pratt/expressions.js");
const { Lexer } = require("./parsing/pratt/lexer.js");
const { ExpressionParser } = require("./parsing/pratt/expression_parser.js");
const { Complex, complex } = require("./math/complex.js");
const { Matrix, matrix } = require("./math/matrix.js");
const { Euclid, Poincare } = require("./math/geometry.js");
const { parameterizePoints, fourierCoefficients, fourierSeries } = require("./math/fourier.js");
const { sscale, ssub, icosphere } = require("./math/icosphere.js"); // fix this garbage
const { stereographic, inverseStereoProject, perspectiveProject } = require("./math/projection.js");
const { rvec } = require("./math/rvector.js");
const { scope, defaultValueScope, valueScope } = require("./scope.js");
const { evaluate } = require("./evaluator.js");

/** -------------------------------------------------------------- */




p5.disableFriendlyErrors = true; // ridiculous that this is on by default
window.plot = undefined;
window.lastMouseX = undefined;
window.lastMouseY = undefined;


/** MathQuill handling */


const MQ = MathQuill.getInterface(2);
const opsString = Object.keys(scope.builtin).filter(x => x.length > 1).join(" ");
const fields = {};

const menuHTML = (id, error=null) => {
    let imageSrc, displayText;
    if (error === null) {
        imageSrc = "http://localhost:8000/data/settings_transparent.png";
        displayText = "Settings";
    } else {
        imageSrc = "http://localhost:8000/data/error_transparent.png";
        displayText = error;
    }
    return `<div>${id}</div>
    <div style="display:flex;">
        <img src="${imageSrc}" style="width:25px;height:25px;" onclick="displayOverlayMenu(${id});" title="${displayText}"></img>
    </div>`;
};

function displayOverlayMenu(id) {
    const overlay = document.querySelector("#overlay-menu-container");
    overlay.style.display = "block";
    overlay.innerHTML = `
    Settings for expression ${id}
    <hr>
    `;
}

document.addEventListener("mousedown", (event) => {
    const overlay = document.querySelector("#overlay-menu-container");
    if (!overlay.contains(event.target)) {
        overlay.style.display = "none";
    }
});

function addField(parent=null) {
    /** add new math input field. parent: parent element */

    const newDiv = document.createElement("div");
    newDiv.setAttribute("class", "math-input-div");

    const newSpan = document.createElement("span");
    newSpan.setAttribute("class", "math-input");

    const newField = MQ.MathField(newSpan, {});
    newDiv.setAttribute("id", `math-input-div-${newField.id}`);

    const newMenu = document.createElement("div");
    newMenu.setAttribute("class", "math-input-side-menu");
    newMenu.setAttribute("id", `math-input-side-menu-${newField.id}`);
    newMenu.innerHTML = menuHTML(newField.id);
    newDiv.appendChild(newMenu);
    newDiv.appendChild(newSpan);

    if (parent === null) {
        const container = document.querySelector("#math-input-container");
        container.appendChild(newDiv);

        fields[newField.id] = {
            id: newField.id,
            field: newField,
            last: null,
            next: null,
            container: newDiv,
        };
    } else {
        const lastDiv = document.querySelector(`#math-input-div-${parent.id}`);
        lastDiv.after(newDiv);

        fields[newField.id] = {
            id: newField.id,
            field: newField,
            last: parent.field,
            next: parent.next,
            container: newDiv,
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
    delete fields[id];
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


function debounceWrapper(func, interval) {
    let timer;
    return function() {
        const context = this;
        const args = arguments;
        clearTimeout(timer);
        timer = setTimeout(() => func.apply(context, args), interval);
    };
}

function getCallbacks(id) {
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

function fieldEditHandler(mathField) {
    /**
     * potentially for the future, instead of debouncing:
     * send the current expression only to the tokenizer,
     * and if it kinda seems ok (well-formed latex 
     * e.g. no \frac{1}{} sorta stuff) then do a full recalc
     */

    scope.userGlobal = {};
    for (const key in valueScope) {
        if (valueScope.hasOwnProperty(key)) {
            delete valueScope[key];
        }
    }
    Object.assign(valueScope, defaultValueScope);

    for (const id of Object.keys(fields)) {
        const latex = fields[id].field.latex();
        fields[id].latex = latex;
        fields[id].expr = cleanLatex(latex);
    }

    const lexer = new Lexer(null, true);
    lexer.setScope(scope);

    const newIdents = [];
    for (const id of Object.keys(fields)) {
        if (!fields[id].expr.includes("=")) continue;
        
        const assignment = fields[id].expr;
        
        const callbacks = getCallbacks(id);
        tracker.setCallback(callbacks["callback"]);
        tracker.setSuccessCallback(callbacks["successCallback"]);
        tracker.clear();

        lexer.setText(assignment);
        lexer.tokenize();

        const tokens = lexer.getTokens();
        const parser = new ExpressionParser(tokens);
        const result = parser.parseExpression();
        
        if (result !== undefined) {
            const left = result.mLeft;
            const isFunction = left instanceof CallExpression;
            const ident = (isFunction) ? left.mFunction.mName : left.mName;
            if (scope.builtin[ident] !== undefined) {
                tracker.error(`cannot overwrite builtin identifier ${ident}`);
                continue;
            } else if (newIdents.includes(ident)) {
                tracker.error(`multiple definitions for ${ident}`);
                continue;
            } else {
                scope.userGlobal[ident] = {
                    isFunction: isFunction,
                };

                if (isFunction) {
                    // define local variables for the function
                    scope.userGlobal[ident].args = {};
                    for (let arg of left.mArgs) {
                        if (arg instanceof OperatorExpression) {
                            if (arg.mLeft instanceof NameExpression && arg.mOperator === TokenType.COLON) {
                                if (arg.mRight instanceof NameExpression && !!scope.builtin[arg.mRight.mName]?.isType) {
                                    scope.userGlobal[ident].args[arg.mLeft.mName] = arg.mRight.mName;
                                } else {
                                    tracker.error("Invalid type annotation");
                                    continue;
                                }
                            } else {
                                tracker.error("Invalid arguments");
                                continue;
                            }
                        } else if (arg instanceof NameExpression) {
                            // this is where type inference should go, if I decide to implement it.
                            // for now, we assume unspecified arguments are complex
                            scope.userGlobal[ident].args[arg.mName] = "complex";
                        } else {
                            tracker.error("Invalid arguments");
                            continue;
                        }
                    }
                }

                newIdents.push(ident);
            }
        }
    }

    if (!tracker.hasError) {
        lexer.setScope(scope); // the scope may have changed, so it doesn't hurt to do this explicitly
        lexer.setAllowUnboundIdentifiers(false);
        for (const id of Object.keys(fields)) {
            const callbacks = getCallbacks(id);
            tracker.setCallback(callbacks["callback"]);
            tracker.setSuccessCallback(callbacks["successCallback"]);

            const expr = fields[id].expr;
            lexer.setText(expr);
            lexer.tokenize();
            const tokens = lexer.getTokens();

            if (tracker.hasError) {
                fields[id].ast = null;
                continue;
            }

            const parser = new ExpressionParser(tokens);
            const result = parser.parseExpression();

            if (tracker.hasError) {
                fields[id].ast = null;
                continue;
            }

            fields[id].ast = result;
            tracker.clear();
        }
    }

    for (const id of Object.keys(fields)) {
        const ast = fields[id].ast;
        if (ast === null) continue;
        evaluate(ast);
    }

    if (valueScope["f"] !== undefined && !tracker.hasError) {
        plot.clear();
        plot.addPlottable(new DomainColoring(
            (z) => valueScope["f"].call({z:z}),
        ));
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
        enter: (mathField) => { advance(mathField.id, 1); },
        downOutOf: (mathField) => { advance(mathField.id, 1); },
        upOutOf: (mathField) => { advance(mathField.id, -1); },
        deleteOutOf: (direction, mathField) => { if (direction === MQ.L) deleteField(mathField.id); },
        edit: debounceWrapper(fieldEditHandler, 500),
    },
});


/** ******************************** */



function linspace(min, max, n) {
	/* Returns n equally spaced values between min and max (including endpoints) */
	const result = [];
	const range = max - min;
	for (let i=0; i<n; i++) {
		result.push(min + range * i / (n-1));
	}
	return result;
}

function roundTo(x, places) {
	/* Rounds x to the specified number of places after the decimal */
	return Math.round(x * Math.pow(10, places)) / Math.pow(10, places);
}

function randInt(a, b) {
    /*
    Generate random integer between a and b
    */
   return Math.floor(a + (b - a + 1) * 0.9999 *Math.random());
}



class Plot {

    static modes = {
        PLANE: 1,
        SPHERE: 2,
        CUBE: 3,
    };

    constructor(displayWidth, displayHeight, bounds=null, mode=null, displayWindowInfo) {
        this.gridlineSpacing = 1;
        this.boundsChangedSinceLastDraw = false;
        this.displayWindowInfo = this.displayWindowInfo;
        this.configureWindow(displayWidth, displayHeight, bounds, bounds === null);
        this.plottables = [];
        this.polygons = [];
        this.needsUpdate = true;

        this.mode = (mode === null) ? Plot.modes.PLANE : mode;
        this.camera = {
            alpha: 1,
            beta: 0,
            pitch: .786,
            roll: 0,
            yaw: .672,
        };
        this.calculateRotationMatrix();
    }

    configureWindow(newWidth=null, newHeight=null, bounds=null) {
        let widthRatio, heightRatio;
        if (newWidth === null) {
            widthRatio = 1;
            heightRatio = 1;
            newWidth = this.dimensions.re;
            newHeight = this.dimensions.im;
        } else {
            if (this.dimensions === undefined) {
                widthRatio = 1;
                heightRatio = 1;
            } else {
                widthRatio = this.dimensions.re / newWidth;
                heightRatio = this.dimensions.im / newHeight;
            }
        }
        
        this.dimensions = complex(newWidth, newHeight)
        this.aspect = newHeight / newWidth;
        this.halfDimensions = this.dimensions.scale(0.5);

        let fitToSquare = false;
        if (bounds === null) {
            if (this.bounds === undefined) {
                bounds = {
                    xMin: -4,
                    xMax: 4,
                    yMin: -4,
                    yMax: 4
                };
                fitToSquare = true;
            } else {
                // retain ratio
                const centerX = (this.bounds.xMin + this.bounds.xMax) / 2;
                const centerY = (this.bounds.yMin + this.bounds.yMax) / 2;
                const halfRangeX = this.units.re * heightRatio / 2;
                const halfRangeY = this.units.im * widthRatio / 2;
                bounds = {
                    xMin: centerX - halfRangeX,
                    xMax: centerX + halfRangeX,
                    yMin: centerY - halfRangeY,
                    yMax: centerY + halfRangeY,
                }
            }
        }

        this.bounds = bounds;
        this.offset = complex(
            ( this.bounds.xMax + this.bounds.xMin ) / 2,
            ( this.bounds.yMax + this.bounds.yMin ) / 2
        );
        this.units = complex(
            this.bounds.xMax - this.bounds.xMin,
            this.bounds.yMax - this.bounds.yMin
        );
        this.pixelsPerUnit = complex(
            this.dimensions.re / this.units.re,
            this.dimensions.im / this.units.im
        );
        this.gridlineCount = complex(
            this.units.re / this.gridlineSpacing,
            this.units.im / this.gridlineSpacing
        );

        if (fitToSquare) this.fitBoundsToSquare();

        valueScope["realBounds"] = complex(this.bounds.xMin, this.bounds.xMax);
        valueScope["imagBounds"] = complex(this.bounds.yMin, this.bounds.yMax);

        this.needsUpdate = true;
        this.boundsChangedSinceLastDraw = true;

        if (this.displayWindowInfo) {
            console.log("Window configuration");
            console.log("window bounds", this.bounds);
            console.log("offset", this.offset);
            console.log("window size", this.dimensions);
            console.log("gridline count", this.gridlineCount);
        }
    }

    fitBoundsToSquare() {
        const centerY = (this.bounds.yMin + this.bounds.yMax) / 2;
        const halfRangeY = (this.units.re * this.aspect) / 2;
        this.configureWindow(this.dimensions.re, this.dimensions.im, {
            xMin: this.bounds.xMin,
            xMax: this.bounds.xMax,
            yMin: centerY - halfRangeY,
            yMax: centerY + halfRangeY
        });
    }

    setCamera(camera) {
        this.camera.alpha = (camera.alpha === undefined) ? this.camera.alpha : camera.alpha;
        this.camera.beta = (camera.beta === undefined) ? this.camera.beta : camera.beta;
        this.camera.pitch = (camera.pitch === undefined) ? this.camera.pitch : camera.pitch;
        this.camera.yaw = (camera.yaw === undefined) ? this.camera.yaw : camera.yaw;
        this.camera.roll = (camera.roll === undefined) ? this.camera.roll : camera.roll;
        this.needsUpdate = true;
    }

    pan(offset) {
        if (this.mode === Plot.modes.PLANE) {
            this.configureWindow(null, null, {
                xMin: this.bounds.xMin + offset.re,
                xMax: this.bounds.xMax + offset.re,
                yMin: this.bounds.yMin + offset.im,
                yMax: this.bounds.yMax + offset.im
            });
        } else {
            this.setCamera({
                pitch: Math.max(Math.min(this.camera.pitch + offset.im, 0.5 * Math.PI), -0.5 * Math.PI),
                yaw: Math.max(Math.min(this.camera.yaw - offset.re, Math.PI), -Math.PI),
            });
        }
    }

    zoom(factor) {
        const newHalfUnits = this.units.scale(factor / 2);
        const center = complex(
            (this.bounds.xMin + this.bounds.xMax) / 2,
            (this.bounds.yMin + this.bounds.yMax) / 2,
        );

        this.configureWindow(null, null, {
            xMin: center.re - newHalfUnits.re,
            xMax: center.re + newHalfUnits.re,
            yMin: center.im - newHalfUnits.im,
            yMax: center.im + newHalfUnits.im,
        });
    }

    state() {
        const latex = [];
        for (const id of Object.keys(fields)) {
            if (fields[id].latex) latex.push(fields[id].latex);
        }
        return JSON.stringify({
            camera: this.camera,
            bounds: this.bounds,
            expressions: latex,
            mode: this.mode,
        });
    }

    loadState(state) {
        state = JSON.parse(state);

        this.setCamera(state.camera);
        this.configureWindow(null, null, state.bounds);
        this.setMode(state.mode);

        for (const id of Object.keys(fields)) {
            deleteField(id, false);
        }

        let lastField = null;
        for (const expr of state.expressions) {
            const newField = addField(lastField);
            fields[newField].field.latex(expr);
        }
    }

    unitsToPixels(z) {
        return complex(
            (z.re - this.offset.re) * this.pixelsPerUnit.re + this.halfDimensions.re,
            -(z.im - this.offset.im) * this.pixelsPerUnit.im + this.halfDimensions.im
        );
    }

    applyCamera(z) {
        if (this.mode === Plot.modes.SPHERE) {
            if (z instanceof Complex) {
                z = inverseStereoProject(z);
            } else {
                z = matrix(z).transpose();
            }
            return Matrix.multiply(this.rotationMatrix, z);
        } else if (this.mode === Plot.modes.CUBE) {
            if (z instanceof Complex) {
                z = matrix([
                    z.re, z.im, 0,
                ]).transpose();
            } else {
                z = matrix(z).transpose();
            }
            return Matrix.multiply(this.rotationMatrix, z);
        } else {
            return z;
        }
    }

    coordinateTransform(z) {
        if (this.mode !== Plot.modes.PLANE) z = perspectiveProject(z, this.camera.alpha, this.camera.beta);
        return this.unitsToPixels(z);
    }

    spaceToScreen(z) {
        return this.coordinateTransform(this.applyCamera(z));
    }

    pixelsToUnits(z) {
        // note: this is NOT an exact inverse of unitsToPixels!!!
        return complex(
            z.re / this.pixelsPerUnit.re,
            -z.im / this.pixelsPerUnit.im
        );
    }

    setMode(mode) {
        this.mode = mode;

        if (mode === Plot.modes.PLANE) {
            if (this.savedBounds !== undefined) this.configureWindow(null, null, this.savedBounds);
        } else {
            this.savedBounds = this.bounds;
            this.configureWindow(null, null, {
                xMin: -2.5,
                xMax: 2.5,
                yMin: -1.5,
                yMax: 1.5,
            });
            this.fitBoundsToSquare();
        }

        this.boundsChangedSinceLastDraw = true;
        this.needsUpdate = true;
    }

    calculateRotationMatrix() {
        this.rotationMatrix = Matrix.rotationMatrix3D(this.camera.pitch, this.camera.roll, this.camera.yaw);
    }

    clear() {
        this.plottables = [];
        this.needsUpdate = true;
    }

    addPlottable(plottable) {
        this.plottables.push(plottable);
        this.needsUpdate = true;
    }

    drawAxes() {
        if (this.mode === Plot.modes.PLANE) {
            const xAxisStart = this.spaceToScreen(complex(this.bounds.xMin, 0));
            const xAxisStop = this.spaceToScreen(complex(this.bounds.xMax, 0));
            const yAxisStart = this.spaceToScreen(complex(0, this.bounds.yMin));
            const yAxisStop = this.spaceToScreen(complex(0, this.bounds.yMax));

            push();
            
            stroke(0);
            strokeWeight(1);
            line(xAxisStart.re, xAxisStart.im, xAxisStop.re, xAxisStop.im);
            line(yAxisStart.re, yAxisStart.im, yAxisStop.re, yAxisStop.im);

            pop();
        } else {
            return;
            const xAxisStart = this.spaceToScreen([-2, 0, 0]);
            const xAxisStop = this.spaceToScreen([2, 0, 0]);
            const yAxisStart = this.spaceToScreen([0, -2, 0]);
            const yAxisStop = this.spaceToScreen([0, 2, 0]);
            const zAxisStart = this.spaceToScreen([0, 0, -2]);
            const zAxisStop = this.spaceToScreen([0, 0, 2]);
            
            push();
            
            stroke(0);
            strokeWeight(1);
            line(xAxisStart.re, xAxisStart.im, xAxisStop.re, xAxisStop.im);
            line(yAxisStart.re, yAxisStart.im, yAxisStop.re, yAxisStop.im);
            line(zAxisStart.re, zAxisStart.im, zAxisStop.re, zAxisStop.im);

            pop();
        }
    }

    drawGridlines() {
        const minBound = (x) => floor(x/this.gridlineSpacing) * this.gridlineSpacing;
        const minBoundX = minBound(this.bounds.xMin);
        const minBoundY = minBound(this.bounds.yMin);

        push();

        stroke(200);
        strokeWeight(1);
        for (let i=0; i<this.gridlineCount.im+1; i++) {
            // horizontal gridlines
            const y = minBoundY + i * this.gridlineSpacing;
            const start = this.unitsToPixels(complex(this.bounds.xMin, y));
            const end = this.unitsToPixels(complex(this.bounds.xMax, y));
            line(start.re, start.im, end.re, end.im);
        }
        for (let i=0; i<this.gridlineCount.re+1; i++) {
            // vertical gridlines
            const x = minBoundX + i * this.gridlineSpacing;
            const start = this.unitsToPixels(complex(x, this.bounds.yMin));
            const end = this.unitsToPixels(complex(x, this.bounds.yMax));
            line(start.re, start.im, end.re, end.im);
        }

        pop();
    }

    draw() {
        background(255);
        this.drawGridlines();
        this.drawAxes();

        this.polygons = [];

        for (let plottable of this.plottables) {
            plottable.update();
            this.polygons = this.polygons.concat(plottable.getPolygons());
        }

        if (this.mode !== Plot.modes.PLANE) {
            this.polygons.sort((poly1, poly2) => {
                return this.applyCamera(poly2.centroid).get(1, 0) - this.applyCamera(poly1.centroid).get(1, 0);
            });
        }

        push();
        for (let poly of this.polygons) {
            if (poly.vertices.length === 1) {
                // point
                fill(0);
                noStroke();
                const point = this.coordinateTransform(this.applyCamera(poly.vertices[0]));
                circle(point.re, point.im, 1);
            } else if (poly.vertices.length === 2) {
                // line
                stroke(0);
                const point1 = this.coordinateTransform(this.applyCamera(poly.vertices[0]));
                const point2 = this.coordinateTransform(this.applyCamera(poly.vertices[1]));
                line(point1.re, point1.im, point2.re, point2.im);
            } else {
                // polygon
                if (poly.outline) {
                    stroke(0);
                } else {
                    noStroke();
                }
                fill(poly.fillColor);

                beginShape();
                for (let vert of poly.vertices) {
                    const point = this.coordinateTransform(this.applyCamera(vert));
                    vertex(point.re, point.im);
                }
                endShape(CLOSE);
            }
        }
        pop();

        this.boundsChangedSinceLastDraw = false;
    }

    update() {
        if (this.needsUpdate) {
            if (this.mode !== Plot.modes.PLANE) this.calculateRotationMatrix();
            this.draw();
            this.needsUpdate = false;
        }
    }

}


class Plottable {

    static id = 0;

    constructor() {
        this.id = Plottable.id;
        Plottable.id++;
    }

    getPolygons() {

    }

    update() {

    }

}


class Point extends Plottable {

    constructor(position) {
        super();
        this.position = position;
        this.radius = 1;
    }

    getPolygons() {
        return [new Polygon([this.position])];
    }

}


class Circle extends Plottable {

    constructor(center, radius) {
        super();
        this.center = center;
        this.radius = radius;
        this.generatePoints();
    }

    generatePoints() {
        this.points = [];
        const pointCount = 100;
        for (let i=0; i<pointCount; i++) {
            const angle = 2 * Math.PI * i / pointCount;
            this.points.push(
                complex(
                    Math.cos(angle), Math.sin(angle)
                ).scale(this.radius).add(this.center)
            );            
        }
    }

    getPolygons() {
        const polygons = [];
        for (let i=1; i<this.points.length; i++) {
            polygons.push(new Polygon(
                [this.points[i-1], this.points[i]]
            ));
        }
        return polygons;
    }

}


class Parametric extends Plottable {

    constructor(fn, range=null, pointCount=null) {
        /**
         * fn: parameterization U->C of curve in C, U subs R
         * U = [range.start, range.stop]
         */
        super();
        this.fn = fn;
        this.pointCount = (pointCount === null) ? 100 : pointCount;
        this.setRange(range);
    }

    setRange(range) {
        this.range = (range === null) ? {start: 0, stop: 1} : range;
        this.generatePoints();
    }

    generatePoints() {
        this.points = [];
        for (let t of linspace(this.range.start, this.range.stop, this.pointCount)) {
            this.points.push(
                this.fn(t)
            );
        }
    }

    getPolygons() {
        const polygons = [];
        for (let i=1; i<this.points.length; i++) {
            polygons.push(new Polygon(
                [this.points[i-1], this.points[i]]
            ));
        }
        return polygons;
    }

}


class Polygon {

    constructor(vertices, fillColor, outline=false) {
        this.vertices = vertices;
        this.fillColor = fillColor;
        this.outline = outline;

        if (vertices[0] instanceof Complex) {
            this.centroid = Euclid.centroid(vertices);
        } else {
            let totalX = 0, totalY = 0, totalZ = 0;
            for (let vert of vertices) {
                totalX += vert[0];
                totalY += vert[1];
                totalZ += vert[2];
            }
            this.centroid = [
                totalX / this.vertices.length,
                totalY / this.vertices.length,
                totalZ / this.vertices.length,
            ];
        }
    }

}


class NormPlot extends Plottable {

    constructor(fn, bounds=null, density=100) {
        super();
        this.fn = fn;
        if (bounds === null) {
            this.bounds = plot.bounds
            this.fixedBounds = false;
        } else {
            this.bounds = bounds;
            this.fixedBounds = true;
        }
        this.samples = complex(density, density);
        this.generatePolygons();
    }

    generatePolygons() {
        this.polygons = [];

        if (plot.mode !== Plot.modes.CUBE) {
            return;
        }

        let x = this.bounds.xMin, y = this.bounds.yMin;
        const step = complex(
            (this.bounds.xMax - this.bounds.xMin) / (this.samples.re - 1),
            (this.bounds.yMax - this.bounds.yMin) / (this.samples.im - 1),
        );

        push(); // for colormode
        colorMode(RGB);
        for (let i=0; i<this.samples.re-1; i++) {
            for (let j=0; j<this.samples.im-1; j++) {
                const square = [
                    complex(x, y),
                    complex(x + step.re + 0.01, y),
                    complex(x + step.re + 0.01, y + step.im + 0.01),
                    complex(x, y + step.im + 0.01),
                ].map(z => [
                    z.re, z.im, this.fn(z).norm(),
                ]);
                const centroid = complex(x + step.re / 2, y + step.im / 2);
                const output = this.fn(centroid);
                const color1 = color(100, 0, 255*Math.tanh(  ((2 * Math.PI + output.arg()) % (2*Math.PI)) / (2 * Math.PI)  ));

                this.polygons.push(new Polygon(square, color1));

                x += step.re;
            }
            x = this.bounds.xMin;
            y += step.im;
        }
        pop();
    }

    getPolygons() {
        return this.polygons;
    }

    update() {
        if (plot.boundsChangedSinceLastDraw && !this.fixedBounds) {
            // TODO: This can be optimized to only recalculate the new polygons
            // by checking the difference between this.bounds and plot.bounds
            this.bounds = plot.bounds;
            this.generatePolygons();
        }
    }

}


class DomainColoring extends Plottable {

    constructor(fn, bounds=null, density=100) {
        super();
        this.fn = fn;
        if (bounds === null) {
            this.bounds = plot.bounds
            this.fixedBounds = false;
        } else {
            this.bounds = bounds;
            this.fixedBounds = true;
        }
        this.samples = complex(density, density);
        this.subdivisions = Math.floor(Math.log(density * density) / Math.log(4));
        this.generatePolygons();
    }

    generatePolygons() {
        if (plot.mode === Plot.modes.PLANE) {
            this.generatePolygonsPlane();
        } else if (plot.mode === Plot.modes.SPHERE) {
            this.generatePolygonsSphere();
        } else {
            this.generatePolygonsPlane();
        }
    }

    generatePolygonsSphere() {
        this.polygons = icosphere(4);
        const threshold = 1000;
        const angleTransform = (angle) => {
            return 360 * ((angle + 2 * Math.PI) % (2 * Math.PI)) / (2 * Math.PI); // modulo is not true remainder in JS
        };
        const normTransform = (norm) => {
            return 25 + 75 * (
                (2 / Math.PI) * Math.atan(norm)
            );
        };
        const parabolaStep = (x) => {
            return x * x * x * (10 - 15 * x + 6 * x * x);
        };
        const highlightPoles = (norm) => {
            return 100 - 85 * Math.max(0, Math.min(1, 0.5 * (Math.sign(norm - threshold) + 1) * parabolaStep(norm - threshold)));
        };

        const filterNaN = (z) => {
            return complex(
                (z.re < z.re + 1) ? z.re : 0,
                (z.im < z.im + 1) ? z.im : 0,
            );
        };

        const getColor = (z) => {
            const zNormal = filterNaN(complex(z.re - Math.floor(z.re), 1 - (z.im - Math.floor(z.im))).eMult(complex(cImage.width-1, cImage.height-1)));
            return cImage.get(zNormal.re, zNormal.im);
        };

        push();
        colorMode(HSB);
        for (let i=0; i<this.polygons.length; i++) {
            this.polygons[i] = new Polygon(this.polygons[i]);
            const centroid = this.polygons[i].centroid;

            // scale slightly from center
            this.polygons[i].vertices = [
                ssub(1, sscale(1.1, ssub(-1, this.polygons[i].vertices[0], centroid)), centroid),
                ssub(1, sscale(1.1, ssub(-1, this.polygons[i].vertices[1], centroid)), centroid),
                ssub(1, sscale(1.1, ssub(-1, this.polygons[i].vertices[2], centroid)), centroid),
            ];

            const output = this.fn(stereographic(centroid));
            const norm = output.norm();
            if (Complex.infinite(output) || Complex.nan(output)) {
                this.polygons[i].fillColor = color(0, 0, 100);
            } else {
                this.polygons[i].fillColor = color(angleTransform(output.arg()), highlightPoles(norm), normTransform(norm));
                // this.polygons[i].fillColor = getColor(output);
            }
        }
        pop();
    }

    generatePolygonsPlane() {
        let x = this.bounds.xMin, y = this.bounds.yMin;
        const step = complex(
            (this.bounds.xMax - this.bounds.xMin) / (this.samples.re - 1),
            (this.bounds.yMax - this.bounds.yMin) / (this.samples.im - 1),
        );
        const angleTransform = (angle) => {
            return 360 * ((angle + 2 * Math.PI) % (2 * Math.PI)) / (2 * Math.PI); // modulo is not true remainder in JS
        };
        const normTransform = (norm) => {
            return 25 + 75 * (
                Math.floor((2 / Math.PI * Math.atan(Math.sqrt(norm))) / 0.2) * 0.2
            );
        };
        const filterNaN = (z) => {
            return complex(
                (z.re < z.re + 1) ? z.re : 0,
                (z.im < z.im + 1) ? z.im : 0,
            );
        };

        const a = 2;
        const distFromAx = (z) => {
            z = rvec(z.re, z.im);
            return Math.tanh(z.proj(rvec(1, a)).mag()) * 2 * Math.PI;
        };
        const angleize = (z) => {
            return Math.tanh(z.norm()) * 2 * Math.PI;
        }

        const getColor = (z) => {
            const zNormal = filterNaN(complex(z.re - Math.floor(z.re), 1 - (z.im - Math.floor(z.im))).eMult(complex(cImage.width-1, cImage.height-1)));
            return cImage.get(zNormal.re, zNormal.im);
        };
        this.polygons = [];
        push();
        colorMode(HSB);
        for (let i=0; i<this.samples.re-1; i++) {
            for (let j=0; j<this.samples.im-1; j++) {
                const square = [
                    complex(x, y),
                    complex(x + step.re + 0.01, y),
                    complex(x + step.re + 0.01, y + step.im + 0.01),
                    complex(x, y + step.im + 0.01),
                ];
                const centroid = complex(x + step.re / 2, y + step.im / 2);
                const output = this.fn(centroid);
                let color1;
                if (Complex.infinite(output) || Complex.nan(output)) {
                    color1 = color(0, 0, 100);                    
                } else {
                    color1 = color(angleTransform(output.arg()), 100, normTransform(output.norm()));
                    // color1 = getColor(output);

                    // const aDist = distFromAx(output);
                    // color1 = color(angleTransform(aDist), 100, 100);
                    // color1 = color(angleTransform(angleize(output)), 100, 100);
                }

                this.polygons.push(new Polygon(square, color1));

                x += step.re;
            }
            x = this.bounds.xMin;
            y += step.im;
        }
        pop();
    }

    getPolygons() {
        return this.polygons;
    }

    update() {
        if (plot.boundsChangedSinceLastDraw && !this.fixedBounds) {
            // TODO: This can be optimized to only recalculate the new polygons
            // by checking the difference between this.bounds and plot.bounds
            this.bounds = plot.bounds;
            this.generatePolygons();
        }
    }

}


class Model extends Plottable {

    constructor(triangles) {
        super();
        this.triangles = [];
        for (let triangle of triangles) {
            this.triangles.push(new Polygon(triangle, color(255, 255, 255), true));
        }
    }

    getPolygons() {
        return this.triangles;
    }

}


function wheelHandler(event) {
    event.preventDefault();
    const factor = 1 + Math.tanh(event.deltaY / 100) / 4;
    plot.zoom(factor);
}


let cImage;
function preload() {
    cImage = loadImage("http://localhost:8000/data/grid_3.png");
}

function setup() {
    const canvasDiv = document.querySelector("#canvas-div");
	const canvas = createCanvas(canvasDiv.offsetWidth, canvasDiv.offsetHeight);
	canvas.parent("canvas-div");
    document.querySelector("#canvas-div").onwheel = wheelHandler;
    plot = new Plot(width, height, null, Plot.modes.SPHERE, false);
    tabSwitch(plot.mode-1);
    // const circ = new Parametric(
    //     t => complex(Math.cos(t), Math.sin(t)),
    //     {start: 0, stop: 2 * Math.PI},
    // );
    // plot.addPlottable(circ);

    /** Domain coloring example */
    // const f = (z) => {
    //     return z.mobius(
    //         complex(1, 0),
    //         complex(0, -1),
    //         complex(1, 0),
    //         complex(0, 1),
    //     );
    // };
    // const f = (z) => {
    //     // return Complex.exp(z);
    //     // return z;
    //     // return Complex.sqrt(z);
    //     // return Complex.mult(z, z);
    //     // return Complex.pow(z, complex(5, 0)).sub(complex(1, 0));
    //     return Complex.cos(z);
    //     // return Complex.sqrt(z.mult(z).scale(-4).sub(complex(1, 0))).sub(Complex.mult(complex(0, 2), z));
    // };
    // const dcPlot = new DomainColoring(f);
    // plot.addPlottable(dcPlot);

    /** Example: chaos game */
    // const maxPoints = 10000;
    // let point = complex(0, 0);
    // const vertices = [];
    // const p = 5;
    // const rad = Poincare.regPolyDist(p, 100);
    // for (let j=0; j<p; j++) {
    //     vertices.push(complex(0, j / p * 2 * Math.PI).exp().scale(rad));
    // }
    // const phi = 2 / (1 + Math.sqrt(5));
    // for (let i=0; i<maxPoints; i++) {
    //     point = Euclid.lerp(point, vertices[randInt(0, p-1)], phi);
    //     // point = Poincare.segment(phi, point, vertices[randInt(0, p-1)]);

    //     plot.addPlottable(new Point(
    //         point
    //     ));
    // }
    // vertices.push(vertices[0]);
    // plot.addPlottable(new Parametric(parameterizePoints(vertices)));

    /** Fourier Series Example */
    // fetch("./data/points.json")
    //     .then((response) => response.json())
    //     .then((json) => {
    //         const points = json.points1;
    //         const result = [];
    //         for (let i=0; i<points["x"].length; i++) {
    //             result.push(complex(points["x"][i], points["y"][i]));
    //         }
    //         const f = parameterizePoints(result);
    //         const para = new Parametric(
    //             f,
    //             {start: 0, stop: 1},
    //         );

    //         const fourierTest = new Parametric(
    //             fourierSeries(f, 4),
    //             {start: 0, stop: 1},
    //             200
    //         );

    //         for (let point of result) {
    //             plot.addPlottable(new Point(point));
    //         }
    //         plot.addPlottable(para);
    //         plot.addPlottable(fourierTest);
    //     });

    /** Example: icosphere model */
    // const icosphereTris = icosphere(3);
    // const sphere = new Model(icosphereTris);
    // plot.addPlottable(sphere);

    lastMouseX = mouseX;
    lastMouseY = mouseY;
}

function mouseDragged() {
	if ((0 <= mouseX && mouseX <= width) && (0 <= mouseY && mouseY <= height)) {
        plot.pan(
            plot.pixelsToUnits(complex(
                lastMouseX - mouseX,
                lastMouseY - mouseY
            ))
        );
		lastMouseX = mouseX;
		lastMouseY = mouseY;
	}
}

function mousePressed() {
	lastMouseX = mouseX;
	lastMouseY = mouseY;
}

function mouseReleased() {
	lastMouseX = 0;
	lastMouseY = 0;
}

function windowResized() {
    const canvasDiv = document.querySelector("#canvas-div");
    resizeCanvas(window.innerWidth * 0.75, document.querySelector("#ui-container").offsetHeight);
    plot.configureWindow(width, height);
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

function draw() {
    plot.update();
}


// there might be a better way to do this, but it's actually fine
window.preload = preload;
window.setup = setup;
window.mouseDragged = mouseDragged;
window.mousePressed = mousePressed;
window.mouseReleased = mouseReleased;
window.windowResized = windowResized;
window.tabSwitch = tabSwitch;
window.draw = draw;

},{"./evaluator.js":1,"./math/complex.js":2,"./math/fourier.js":3,"./math/geometry.js":4,"./math/icosphere.js":5,"./math/matrix.js":6,"./math/projection.js":7,"./math/rvector.js":8,"./parsing/errors.js":9,"./parsing/latex_convert.js":10,"./parsing/pratt/expression_parser.js":11,"./parsing/pratt/expressions.js":12,"./parsing/pratt/lexer.js":13,"./parsing/pratt/tokentype.js":18,"./scope.js":20}]},{},[21]);
