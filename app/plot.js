const { Complex, complex } = require("../math/complex.js");
const { Matrix, matrix } = require("../math/matrix.js");
const { icosphere_flat } = require("../math/icosphere.js");
const { inverseStereoProject, perspectiveProject } = require("../math/projection.js");
const { GRADIENTS } = require("./coloring.js");
const { loadImage } = require("./data_loading.js");
const { generateSettingsHTML, addField, deleteField, tabSwitch } = require("./domhandler.js");


const pValueArray = [
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


class Plot {

    static modes = {
        PLANE: 1,
        SPHERE: 2,
        CUBE: 3,
    };

    constructor(
        displayWidth, displayHeight, bounds=null, mode=null, displayWindowInfo, reglInstance=null,
        shaders=null
    ) {
        this.gridlineSpacing = 1;
        this.boundsChangedSinceLastDraw = false;
        this.imageFile = null;
        this.displayWindowInfo = displayWindowInfo;
        this.configureWindow(displayWidth, displayHeight, bounds);
        this.plottables = [];
        this.polygons = [];

        this.reglInstance = reglInstance;
        this.shaders = shaders;
        this.shaderReplacement = null;
        this.activeGradient = null;
        this.setDisplayReplacement(null);

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

        this.planeMesh = this.generatePlaneMesh(100);
        this.sphereMesh = icosphere_flat(5);
    }

    configureWindow(newWidth, newHeight, bounds=null) {
        if (bounds === null) {
            if (!this.windowConfigured) {
                // set to default of 8 real units, centered, and square unit aspect ratio
                this.units = complex(8, 8 * newHeight / newWidth);
                this.offset = complex(0, 0); // centered
            } else {
                // retain ratio
                // oldX / oldWidth = newX / newWidth --> newX = newWidth * (oldX / oldWidth)
                // newY = newHeight * (oldY / oldHeight)
                this.units = complex(
                    newWidth * (this.units.re / this.dimensions.re),
                    newHeight * (this.units.im / this.dimensions.im)
                );
            }
            this.bounds = {
                xMin: this.offset.re - 0.5 * this.units.re,
                xMax: this.offset.re + 0.5 * this.units.re,
                yMin: this.offset.im - 0.5 * this.units.im,
                yMax: this.offset.im + 0.5 * this.units.im,
            };
        } else {
            // set bounds directly
            this.bounds = bounds;
            this.units = complex(this.bounds.xMax - this.bounds.xMin, this.bounds.yMax - this.bounds.yMin);
            this.offset = complex(
                (this.bounds.xMin + this.bounds.xMax) / 2,
                (this.bounds.yMin + this.bounds.yMax) / 2,
            );
        }
        this.dimensions = complex(newWidth, newHeight);
        this.aspect = newHeight / newWidth;
        this.halfDimensions = this.dimensions.scale(0.5);
        this.pixelsPerUnit = complex(
            this.dimensions.re / this.units.re,
            this.dimensions.im / this.units.im
        );
        this.gridlineCount = complex(
            this.units.re / this.gridlineSpacing,
            this.units.im / this.gridlineSpacing
        );

        this.boundsChangedSinceLastDraw = true;
        this.needsUpdate = true;
        this.windowConfigured = true;
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
            this.configureWindow(this.dimensions.re, this.dimensions.im, {
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

        this.configureWindow(this.dimensions.re, this.dimensions.im, {
            xMin: center.re - newHalfUnits.re,
            xMax: center.re + newHalfUnits.re,
            yMin: center.im - newHalfUnits.im,
            yMax: center.im + newHalfUnits.im,
        });
    }

    uploadImage(id) {
        const tempEl = document.createElement("input");
        tempEl.setAttribute("type", "file");
        tempEl.setAttribute("accept", "image/*");
        tempEl.setAttribute("id", "file-selector");
        document.body.appendChild(tempEl);
        tempEl.click();

        tempEl.onchange = () => {
            const selector = document.querySelector("#file-selector");
            const files = selector.files;
            if (files?.length <= 0) return;
            const reader = new FileReader();

            reader.onload = (event) => {
                loadImage(URL.createObjectURL(files[0])).then(result => {
                    this.shaders["sample texture"] = result;
                    fields[id]["imageFile"] = files[0].name;
                    fieldEditHandler(null);
                    fields[id]["settingsHTML"] = generateSettingsHTML(id);
                    displayOverlayMenu(id);
                });
            }

            reader.readAsText(files.item(0));
            selector.parentNode.removeChild(selector);
        }
    }

    state() {
        const latex = [];
        for (const id of Object.keys(fields)) {
            if (fields[id]) latex.push(fields[id].field.latex());
        }
        return JSON.stringify({
            camera: this.camera,
            bounds: this.bounds,
            expressions: latex,
            mode: this.mode,
        }, null, 4);
    }

    loadState(state) {
        state = JSON.parse(state);

        this.setCamera(state.camera);
        this.configureWindow(this.dimensions.re, this.dimensions.im, state.bounds);
        tabSwitch(state.mode-1);

        for (const id of Object.keys(fields)) {
            deleteField(id, false);
        }

        let lastField = null;
        for (const expr of state.expressions) {
            const newField = addField(lastField);
            fields[newField].field.latex(expr);
        }
    }

    downloadState() {
        const state = encodeURIComponent(this.state());
        const tempEl = document.createElement("a");
        tempEl.setAttribute("href", `data:text/json;charset=utf-8,${state}`);
        tempEl.setAttribute("download", "plot.json");
        document.body.appendChild(tempEl);
        tempEl.click();
        tempEl.remove();
    }

    uploadState() {
        const tempEl = document.createElement("input");
        tempEl.setAttribute("type", "file");
        tempEl.setAttribute("accept", ".json");
        tempEl.setAttribute("id", "file-selector");
        document.body.appendChild(tempEl);
        tempEl.click();

        tempEl.onchange = () => {
            const selector = document.querySelector("#file-selector");
            const files = selector.files;
            if (files?.length <= 0) return;
            const reader = new FileReader();

            reader.onload = (event) => {
                this.loadState(event.target.result);
            }

            reader.readAsText(files.item(0));
            selector.parentNode.removeChild(selector);
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
        const previousMode = this.mode;
        this.mode = mode;

        if (mode === Plot.modes.PLANE) {
            if (this.savedBounds !== undefined) this.configureWindow(this.dimensions.re, this.dimensions.im, this.savedBounds);
        } else {
            if (previousMode === Plot.modes.PLANE) this.savedBounds = this.bounds;
            this.configureWindow(this.dimensions.re, this.dimensions.im, {
                xMin: -2.5,
                xMax: 2.5,
                yMin: -2.5 * this.aspect,
                yMax: 2.5 * this.aspect,
            });
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

    generatePlaneMesh(count) {
        let x = -1, y = -1;
        const step = 2 / count;
        const verts = [];
        const sqrt2 = Math.sqrt(2);
        for (let i=0; i<count; i++) {
            for (let j=0; j<count; j++) {
                const nextX = x + step;
                const nextY = y + step;

                verts.push([x / sqrt2, y / sqrt2, 0]);
                verts.push([nextX / sqrt2, y / sqrt2, 0]);
                verts.push([x / sqrt2, nextY / sqrt2, 0]);

                verts.push([nextX / sqrt2, y / sqrt2, 0]);
                verts.push([nextX / sqrt2, nextY / sqrt2, 0]);
                verts.push([x / sqrt2, nextY / sqrt2, 0]);

                x += step;
            }
            x = -1;
            y += step;
        }
        return verts;
    }

    setShaderReplacement(glslSource) {
        this.shaderReplacement = glslSource;
        this.needsUpdate = true;
    }

    setDisplayReplacement(name, colorGLSL) {
        if (name) {
            this.displayReplacement = `vec2 outp = udf_${name}(z);\n${colorGLSL}`;
            this.displayReplacementFunction = `vec2 outp = udf_${name}(z);`
        } else {
            this.displayReplacement = "vec2 outp = vec2(1., 0.);vec3 col=vec3(0.9, 0.9, 0.9);";
            this.displayReplacementFunction = "vec2 outp = vec2(1., 0.);";
        }
        this.needsUpdate = true;
    }

    toggleDisplay(id) {
        const checkbox = document.querySelector(`#display-checkbox-${id}`);
        fields[id]["displaySettings"]["display"] = checkbox.checked;
        fieldEditHandler(null);
    }

    setColorMode(id) {
        const select = document.querySelector(`#display-coloring-dropdown-${id}`);
        fields[id]["displaySettings"]["colorMode"] = select.value;
        fields[id]["settingsHTML"] = generateSettingsHTML(id);
        displayOverlayMenu(id);
        fieldEditHandler(null);
    }

    setGradientMode(id) {
        const select = document.querySelector(`#gradient-dropdown-${id}`);
        fields[id]["displaySettings"]["gradient"] = select.value;
        this.activeGradient = select.value;
        fieldEditHandler(null);
    }

    drawFnPlane() {
        let frag = this.shaders["complexmos.frag"];
        if (this.shaderReplacement !== null) {
            frag = frag.replace(/\/\/REPLACE_BEGIN.*\/\/REPLACE_END/ms, this.shaderReplacement);
        }
        frag = frag.replace(/\/\/DISPLAY_REPLACE_BEGIN.*\/\/DISPLAY_REPLACE_END/ms, this.displayReplacement);
        const activeGradient = (this.activeGradient) ? this.activeGradient : "monokai";

        return this.reglInstance({
            frag: frag,
            vert: this.shaders["complexmos.vert"],

            attributes: {
                position: [
                    [-1, -1], [1, 1], [-1, 1],
                    [-1, -1], [1, 1], [1, -1],
                ],
            },

            uniforms: {
                width: this.reglInstance.context("viewportWidth"),
                height: this.reglInstance.context("viewportHeight"),

                xBounds: [this.bounds.xMin, this.bounds.xMax],
                yBounds: [this.bounds.yMin, this.bounds.yMax],

                pValues: pValueArray,

                gradR: GRADIENTS[activeGradient][0],
                gradG: GRADIENTS[activeGradient][1],
                gradB: GRADIENTS[activeGradient][2],

                texture: this.reglInstance.texture(this.shaders["sample texture"]),
            },

            count: 6
        });
    }

    drawFnSphere() {
        const mesh = this.sphereMesh;
        const vertexCount = mesh.length;

        let frag = this.shaders["complexmos_sphere.frag"];
        if (this.shaderReplacement !== null) {
            frag = frag.replace(/\/\/REPLACE_BEGIN.*\/\/REPLACE_END/ms, this.shaderReplacement);
        }
        frag = frag.replace(/\/\/DISPLAY_REPLACE_BEGIN.*\/\/DISPLAY_REPLACE_END/ms, this.displayReplacement);
        const activeGradient = (this.activeGradient) ? this.activeGradient : "monokai";

        return this.reglInstance({
            frag: frag,
            vert: this.shaders["complexmos_sphere.vert"],

            attributes: {
                position: mesh,
            },

            uniforms: {
                width: this.reglInstance.context("viewportWidth"),
                height: this.reglInstance.context("viewportHeight"),

                xBounds: [this.bounds.xMin, this.bounds.xMax],
                yBounds: [this.bounds.yMin, this.bounds.yMax],

                pValues: pValueArray,
                alpha: this.camera.alpha,
                beta: this.camera.beta,

                row1: this.rotationMatrix.getRow(0),
                row2: this.rotationMatrix.getRow(1),
                row3: this.rotationMatrix.getRow(2),

                gradR: GRADIENTS[activeGradient][0],
                gradG: GRADIENTS[activeGradient][1],
                gradB: GRADIENTS[activeGradient][2],

                texture: this.reglInstance.texture(this.shaders["sample texture"]),
            },

            count: vertexCount,

            cull: {
                enable: true,
                face: "back",
            },

            frontFace: "cw",
        });
    }

    drawFn3D() {
        const mesh = this.planeMesh;
        const vertexCount = mesh.length;

        let frag = this.shaders["complexmos_cube.frag"];
        let vert = this.shaders["complexmos_cube.vert"];
        if (this.shaderReplacement !== null) {
            frag = frag.replace(/\/\/REPLACE_BEGIN.*\/\/REPLACE_END/ms, this.shaderReplacement);
            vert = vert.replace(/\/\/REPLACE_BEGIN.*\/\/REPLACE_END/ms, this.shaderReplacement);
            vert = vert.replace(/\/\/DISPLAY_REPLACE_BEGIN.*\/\/DISPLAY_REPLACE_END/ms, this.displayReplacementFunction);
        }
        frag = frag.replace(/\/\/DISPLAY_REPLACE_BEGIN.*\/\/DISPLAY_REPLACE_END/ms, this.displayReplacement);
        const activeGradient = (this.activeGradient) ? this.activeGradient : "monokai";

        return this.reglInstance({
            frag: frag,
            vert: vert,

            attributes: {
                position: mesh,
            },

            uniforms: {
                width: this.reglInstance.context("viewportWidth"),
                height: this.reglInstance.context("viewportHeight"),

                xBounds: [this.bounds.xMin, this.bounds.xMax],
                yBounds: [this.bounds.yMin, this.bounds.yMax],

                pValues: pValueArray,
                alpha: this.camera.alpha,
                beta: this.camera.beta,

                row1: this.rotationMatrix.getRow(0),
                row2: this.rotationMatrix.getRow(1),
                row3: this.rotationMatrix.getRow(2),

                gradR: GRADIENTS[activeGradient][0],
                gradG: GRADIENTS[activeGradient][1],
                gradB: GRADIENTS[activeGradient][2],

                texture: this.reglInstance.texture(this.shaders["sample texture"]),
            },

            cull: {
                enable: false,
            },

            frontFace: "ccw",

            count: vertexCount,
        });
    }

    update() {
        if (this.needsUpdate) {
            if (this.mode !== Plot.modes.PLANE) this.calculateRotationMatrix();
            if (RENDERER === "p5") {
                this.draw();
            } else {
                let drawFn;
                if (this.mode === Plot.modes.PLANE) {
                    drawFn = this.drawFnPlane();
                } else if (this.mode === Plot.modes.SPHERE) {
                    drawFn = this.drawFnSphere();
                } else {
                    drawFn = this.drawFn3D();
                }
                drawFn();
            }
            this.needsUpdate = false;
        }
    }

}


module.exports = {
    Plot,
};