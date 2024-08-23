const { Euclid } = require("../math/geometry.js");
const { sscale, ssub, icosphere } = require("../math/icosphere.js");
const { stereographic } = require("../math/projection.js");
const { rvec } = require("../math/rvector.js");
const { Complex, complex } = require("../math/complex.js");



function linspace(min, max, n) {
	/* Returns n equally spaced values between min and max (including endpoints) */
	const result = [];
	const range = max - min;
	for (let i=0; i<n; i++) {
		result.push(min + range * i / (n-1));
	}
	return result;
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
            if (output === null || Complex.infinite(output) || Complex.nan(output)) {
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
                if (output === null || Complex.infinite(output) || Complex.nan(output)) {
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

module.exports = {
    Plottable, Point, Circle, Parametric,
    Polygon, NormPlot, DomainColoring, Model,
};