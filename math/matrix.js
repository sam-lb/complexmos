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
		// ahAhA this relies on the fact that the default inspect window in most browsers uses a monospaced font
		let str = this.toString();
		console.log(str);
	}

}

function matrix(arr) {
    // utility function for easy matrix creation
	return new Matrix(arr);
}
