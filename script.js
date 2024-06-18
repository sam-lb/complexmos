p5.disableFriendlyErrors = true; // ridiculous that this is on by default
let plot, lastMouseX, lastMouseY;


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
    valueScope = {};
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
