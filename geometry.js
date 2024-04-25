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