const { complex } = require("../math/complex.js");
const { hideOverlayMenu } = require("./domhandler.js");

window.lastMouseX = undefined;
window.lastMouseY = undefined;
window.mouseIsDown = false;
window.resizeBarStart = false;

function wheelHandler(event) {
    event.preventDefault();
    const factor = 1 + Math.tanh(event.deltaY / 100) / 4;
    plot.zoom(factor);
}

function mouseDragged(event) {
    if (mouseIsDown) {
        const rect = event.target.getBoundingClientRect();
        const mouseX = (event.touches) ? event.touches[0].clientX - rect.left : event.clientX - rect.left;
        const mouseY = (event.touches) ? event.touches[0].clientY - rect.top : event.clientY - rect.top;
        const canvasDiv = document.querySelector("#canvas-div");
        
        if ((0 <= mouseX && mouseX <= canvasDiv.offsetWidth) && (0 <= mouseY && mouseY <= canvasDiv.offsetHeight)) {
            const diff = complex(lastMouseX - mouseX, lastMouseY - mouseY);
            if (plot.mode === Plot.modes.PLANE) {
                plot.pan(plot.pixelsToUnits(diff));
            } else {
                plot.pan(diff.eMult(complex(3 / canvasDiv.clientWidth, -3 / canvasDiv.clientHeight)));
            }
            lastMouseX = mouseX;
            lastMouseY = mouseY;
        }
    }
}

function mousePressed(event) {
    const rect = event.target.getBoundingClientRect();
    const mouseX = (event.touches) ? event.touches[0].clientX - rect.left : event.clientX - rect.left;
    const mouseY = (event.touches) ? event.touches[0].clientY - rect.top : event.clientY - rect.top;
	
    lastMouseX = mouseX;
	lastMouseY = mouseY;
    mouseIsDown = true;
}

function exprBarMousePressed(event) {
    resizeBarStart = true;
}

function mouseReleased(event) {
	lastMouseX = 0;
	lastMouseY = 0;
    mouseIsDown = false;
    resizeBarStart = false;
}

function exprBarResize(event, callback) {
    if (mouseIsDown && resizeBarStart) {
        const target = document.querySelector("#drag-expr-bar");
        const rect = target.getBoundingClientRect();
        const mouseX = (event.touches) ? event.touches[0].clientX - rect.left : event.clientX - rect.left;
        const dx = lastMouseX - mouseX;

        let x = (rect.right - dx) / window.innerWidth;
        x = Math.max(0.20, Math.min(0.75, x));
        const p = 3 * x / (1 - x);

        document.querySelector("#ui-container").style.flex = p.toString();
        callback();
    }
}

function registerMouseEvents() {
    const canvasDiv = document.querySelector("#canvas-div");
    const dragDiv = document.querySelector("#drag-expr-bar");

    canvasDiv.addEventListener("wheel", wheelHandler);
    canvasDiv.addEventListener("mousemove", mouseDragged);
    canvasDiv.addEventListener("touchmove", mouseDragged);
    document.addEventListener("mousedown", mousePressed);
    document.addEventListener("touchstart", mousePressed);
    document.addEventListener("mouseup", mouseReleased);
    document.addEventListener("touchend", mouseReleased);
    dragDiv.addEventListener("mousedown", exprBarMousePressed);
    dragDiv.addEventListener("touchstart", exprBarMousePressed);
    document.addEventListener("mousemove", (event) => exprBarResize(event, resizeDebounced));
    document.addEventListener("touchmove", (event) => exprBarResize(event, resizeDebounced));
}

function windowResized(callback) {
    setTimeout(() => {
        callback(plot.offset);
    }, 100);
}

document.addEventListener("mousedown", (event) => {
    const overlay = document.querySelector("#overlay-menu-container");
    if (!overlay.contains(event.target)) {
        hideOverlayMenu();
    }
});


module.exports = {
    registerMouseEvents, windowResized,
};