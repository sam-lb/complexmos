const EPSILON = 0.000001;


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
		return this.mult(a).add(b).div(this.mult(c).add(d));
	}

	/* ---------- Static functions -------------------- */

	static norm_(z) {
		return z.norm();
	}

	static normSq_(z) {
		return z.normSq();
	}

	static inv_(z) {
		return z.conj().scale(1 / z.normSq());
	}

	static add(z1, z2) {
		return z1.add(z2);
	}

	static mult(z1, z2) {
		return z1.mult(z2);
	}

	static div_(z1, z2) {
		return z1.div(z2);
	}

	static sin_(z) {
		/*
		Return the principle branch of sin(z)
		*/
		return z.sin();
	}

	/* --------------- In-place operations --------------------- */

	iadd(z) {
		this.re += z.re;
		this.im += z.im;
	}


}

function complex(real, imaginary) {
	/* instantiate a Complex without new keyword */
	return new Complex(real, imaginary);
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

	add(Z) {
		const result = [];
		for (let i=0; i<this.dimension; i++) {
			result.push(Complex.add_(this.get(i), Z.get(i)));
		}
		return result;
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

}


function vector(...values) {
	return new Vector(values);
}