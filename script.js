p5.disableFriendlyErrors = true; // ridiculous that this is on by default
let plot, lastMouseX, lastMouseY;



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

    constructor(displayWidth, displayHeight, bounds=null) {
        this.gridlineSpacing = 1;
        this.boundsChangedSinceLastDraw = false;
        this.configureWindow(displayWidth, displayHeight, bounds, bounds === null);
        this.plottables = [];
        this.needsUpdate = true;
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

        this.needsUpdate = true;
        this.boundsChangedSinceLastDraw = true;

        console.log("Window configuration");
        console.log("window bounds", this.bounds);
        console.log("offset", this.offset);
        console.log("window size", this.dimensions);
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

    pan(offset) {
        this.configureWindow(null, null, {
            xMin: this.bounds.xMin + offset.re,
            xMax: this.bounds.xMax + offset.re,
            yMin: this.bounds.yMin + offset.im,
            yMax: this.bounds.yMax + offset.im
        });
    }

    zoom(factor) {
        this.configureWindow(null, null, {
            xMin: this.bounds.xMin * factor,
            xMax: this.bounds.xMax * factor,
            yMin: this.bounds.yMin * factor,
            yMax: this.bounds.yMax * factor,
        });
    }

    unitsToPixels(z) {
        return complex(
            (z.re - this.offset.re) * this.pixelsPerUnit.re + this.halfDimensions.re,
            -(z.im - this.offset.im) * this.pixelsPerUnit.im + this.halfDimensions.im
        );
    }

    pixelsToUnits(z) {
        // note: this is NOT an exact inverse of unitsToPixels!!!
        return complex(
            z.re / this.pixelsPerUnit.re,
            -z.im / this.pixelsPerUnit.im
        );
    }

    addPlottable(plottable) {
        this.plottables.push(plottable);
        this.needsUpdate = true;
    }

    drawAxes() {
        const xAxisStart = this.unitsToPixels(complex(this.bounds.xMin, 0));
        const xAxisStop = this.unitsToPixels(complex(this.bounds.xMax, 0));
        const yAxisStart = this.unitsToPixels(complex(0, this.bounds.yMin));
        const yAxisStop = this.unitsToPixels(complex(0, this.bounds.yMax));

        push();
        
        stroke(0);
        strokeWeight(1);
        line(xAxisStart.re, xAxisStart.im, xAxisStop.re, xAxisStop.im);
        line(yAxisStart.re, yAxisStart.im, yAxisStop.re, yAxisStop.im);

        pop();
    }

    drawGridlines() {
        const minBound = (x) => { return (floor(x / this.gridlineSpacing) - 1) * this.gridlineSpacing; };
        const minBoundX = minBound(this.bounds.xMin);
        const minBoundY = minBound(this.bounds.yMin);

        push();

        stroke(200);
        strokeWeight(1);
        for (let i=0; i<this.gridlineCount.im+2; i++) {
            // horizontal gridlines
            const y = minBoundY + i * this.gridlineSpacing;
            const start = this.unitsToPixels(complex(this.bounds.xMin, y));
            const end = this.unitsToPixels(complex(this.bounds.xMax, y));
            line(start.re, start.im, end.re, end.im);
        }
        for (let i=0; i<this.gridlineCount.im+3; i++) {
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

        for (let plottable of this.plottables) {
            plottable.draw();
        }
        this.boundsChangedSinceLastDraw = false;
    }

    update() {
        if (this.needsUpdate) {
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

    draw() {

    }

}


class Point extends Plottable {

    constructor(position) {
        super();
        this.position = position;
        this.radius = 1;
    }

    draw() {
        push();
        fill(0);
        noStroke();
        const convertedPos = plot.unitsToPixels(this.position);
        circle(convertedPos.re, convertedPos.im, this.radius);
        pop();
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

    draw() {
        push();

        stroke(0);
        strokeWeight(1);
        noFill();
        beginShape();
        for (let point of this.points) {
            const convertedPoint = plot.unitsToPixels(point);
            vertex(convertedPoint.re, convertedPoint.im);
        }
        endShape(CLOSE);

        pop();
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

    draw() {
        push();

        stroke(0);
        strokeWeight(1);
        noFill();
        beginShape();
        for (let point of this.points) {
            const convertedPoint = plot.unitsToPixels(point);
            vertex(convertedPoint.re, convertedPoint.im);
        }
        endShape();

        pop();
    }

}


class Polygon {

    constructor(vertices, fillColor, outline=false) {
        this.vertices = vertices;
        this.fillColor = fillColor;
        this.outline = outline;
    }

    draw() {
        push();

        if (this.outline) {
            stroke(0);
        } else {
            noStroke();
        }
        fill(this.fillColor);

        beginShape();
        for (let vert of this.vertices) {
            vert = plot.unitsToPixels(vert);
            vertex(vert.re, vert.im);
        }
        endShape(CLOSE);

        pop();
    }

}


class DomainColoring extends Plottable {

    constructor(fn, bounds=null) {
        super();
        this.fn = fn;
        if (bounds === null) {
            this.bounds = plot.bounds
            this.fixedBounds = false;
        } else {
            this.bounds = bounds;
            this.fixedBounds = true;
        }
        this.samples = complex(100, 100);
        this.generatePolygons();
    }

    generatePolygons() {
        let x = this.bounds.xMin, y = this.bounds.yMin;
        const step = complex(
            (this.bounds.xMax - this.bounds.xMin) / (this.samples.re - 1),
            (this.bounds.yMax - this.bounds.yMin) / (this.samples.im - 1),
        );
        const angleTransform = (angle) => {
            return 360 * ((angle + 2 * Math.PI) % (2 * Math.PI)) / (2 * Math.PI); // modulo is not true remainder in JS
        }
        const normTransform = (norm) => {
            // return 25 + (150 / Math.PI) * Math.atan(Math.sqrt(norm));
            return 25 + 75 * (
                Math.floor((2 / Math.PI * Math.atan(Math.sqrt(norm))) / 0.2) * 0.2
            );
        }
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
                const color1 = color(angleTransform(output.arg()), 100, normTransform(output.norm()));

                this.polygons.push(new Polygon(square, color1));

                x += step.re;
            }
            x = this.bounds.xMin;
            y += step.im;
        }
        pop();
    }

    draw() {
        if (plot.boundsChangedSinceLastDraw && !this.fixedBounds) {
            // TODO: This can be optimized to only recalculate the new polygons
            // by checking the difference between this.bounds and plot.bounds
            this.bounds = plot.bounds;
            this.generatePolygons();
        }

        for (let poly of this.polygons) {
            poly.draw();
        }
    }

}




function setup() {
    const canvasDiv = document.querySelector("#canvas-div");
	const canvas = createCanvas(canvasDiv.offsetWidth, canvasDiv.offsetHeight);
	canvas.parent("canvas-div");
    plot = new Plot(width, height);
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
    const f = (z) => {
        // return Complex.exp(z);
        // return z;
        // return Complex.sqrt(z);
        // return Complex.mult(z, z);
        // return Complex.pow(z, complex(5, 0)).sub(complex(1, 0));
        return Complex.cos(z);
    };
    const dcPlot = new DomainColoring(f);
    plot.addPlottable(dcPlot);

    /** Example: chaos game */
    // const maxPoints = 100000;
    // let point = complex(0, 0);
    // const vertices = [];
    // const p = 4;
    // const rad = Poincare.regPolyDist(p, 100);
    // for (let j=0; j<p; j++) {
    //     vertices.push(complex(0, j / p * 2 * Math.PI).exp().scale(rad));
    // }
    // const phi = 2 / (1 + Math.sqrt(5));
    // for (let i=0; i<maxPoints; i++) {
    //     // point = Euclid.lerp(point, vertices[randInt(0, p-1)], phi);
    //     point = Poincare.segment(phi, point, vertices[randInt(0, p-1)]);

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
    resizeCanvas(canvasDiv.offsetWidth, canvasDiv.offsetHeight);
    plot.configureWindow(width, height);
}

function draw() {
    plot.update();
}
