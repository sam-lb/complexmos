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
		 * return whether one or both of the components of z is infinite
		 */
		return (
			z.re === Infinity || z.re === -Infinity 
			|| z.im === Infinity || z.im === -Infinity
		);
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