function preload() {
    window.cImage = loadImage("../data/grid_3.png");
}

async function loadImage(src) {
    return new Promise((resolve, reject) => {
        let image = new Image();
        image.src = src;
        image.onerror = reject;
        image.onload = resolve(image);
    });
}

async function loadShaders() {
    const frag = (await fetch("../shaders/complexmos.frag"));
    const vert = (await fetch("../shaders/complexmos.vert"));
    const fragSphere = (await fetch("../shaders/complexmos_sphere.frag"));
    const vertSphere = (await fetch("../shaders/complexmos_sphere.vert"));
    const fragCube = (await fetch("../shaders/complexmos_cube.frag"));
    const vertCube = (await fetch("../shaders/complexmos_cube.vert"));
    const complexLib = (await fetch("../shaders/complex.frag"));
    const coloringLib = (await fetch("../shaders/coloring.frag"));
    const sampleImage = await loadImage("../data/cat.jpg");

    const complexLibSource = await complexLib.text().then(text => text);
    const coloringLibSource = await coloringLib.text().then(text => text);

    const importLib = (match, replacement) => (fileContents) => fileContents.replace(match, replacement);
    const importComplex = importLib(/\/\/IMPORT_COMPLEX/, complexLibSource);
    const importColoring = importLib(/\/\/IMPORT_COLORING/, coloringLibSource);

    const fragShaderSource = await frag.text().then(importComplex).then(importColoring);
    const vertShaderSource = await vert.text().then(text => text);
    const fragCubeShaderSource = await fragCube.text().then(importComplex).then(importColoring);
    const vertCubeShaderSource = await vertCube.text().then(importComplex); // yes the vert shader needs this
    const fragSphereShaderSource = await fragSphere.text().then(importComplex).then(importColoring);
    const vertSphereShaderSource = await vertSphere.text().then(text => text);

    return {
        "complexmos.frag": fragShaderSource,
        "complexmos.vert": vertShaderSource,
        "complexmos_sphere.frag": fragSphereShaderSource,
        "complexmos_sphere.vert": vertSphereShaderSource,
        "complexmos_cube.frag": fragCubeShaderSource,
        "complexmos_cube.vert": vertCubeShaderSource,
        "sample texture": sampleImage,
    };
}

module.exports = {
    loadImage, loadShaders, preload,
};