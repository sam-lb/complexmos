



function snormal(X) {
    const norm = Math.sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
    return [
        X[0] / norm,
        X[1] / norm,
        X[2] / norm,
    ];
}

function ssub(s, X, Y) {
    return [
        X[0] + s * Y[0],
        X[1] + s * Y[1],
        X[2] + s * Y[2],
    ];
}

function sscale(s, X) {
    return [
        s * X[0],
        s * X[1],
        s * X[2],
    ];
}



function icosphere(subdivisions=0) {

    // original icosahedron data

    const scale = 1 / (2 * Math.sin(2 * Math.PI / 5));
    const phi = scale * (1 + Math.sqrt(5)) / 2;

    const vertices = [
        [-scale, phi, 0], [scale, phi, 0], [-scale, -phi, 0], [scale, -phi, 0],
        [0, -scale, phi], [0, scale, phi], [0, -scale, -phi], [0, scale, -phi],
        [phi, 0, -scale], [phi, 0, scale], [-phi, 0, -scale], [-phi, 0, scale],
    ];

    const triangleIndices = [
        [0, 11, 5],[0, 5, 1],[0, 1, 7],[0, 7, 10],[0, 10, 11],
        [11, 10, 2],[5, 11, 4],[1, 5, 9],[7, 1, 8],[10, 7, 6],
        [3, 9, 4],[3, 4, 2],[3, 2, 6],[3, 6, 8],[3, 8, 9],
        [9, 8, 1],[4, 9, 5],[2, 4, 11],[6, 2, 10],[8, 6, 7],
    ];

    let triangles = [];
    for (let index of triangleIndices) {
        triangles.push([vertices[index[0]], vertices[index[1]], vertices[index[2]]]);
    }

    // split into 9:

    // const baryCoefs = [
    //     [0, 0], [1/3, 0], [2/3, 0], [1, 0], [0, 1/3],
    //     [1/3, 1/3], [2/3, 1/3], [0, 2/3], [1/3, 2/3], [0, 1],
    // ];

    // const trindices = [
    //     [0, 1, 4], [1, 2, 5], [2, 3, 6], [4, 5, 1], [5, 6, 2],
    //     [4, 5, 7], [5, 6, 8], [7, 8, 5], [7, 8, 9],
    // ];

    // split into 4:

    const baryCoefs = [
        [0, 0], [.5, 0], [1, 0],
        [0, .5], [.5, .5], [0, 1],
    ];

    const trindices = [
        [0, 1, 3], [1, 2, 4],
        [1, 3, 4], [3, 4, 5],
    ];

    for (let i=0; i<subdivisions; i++) {
        const newTriangles = [];
        for (let tri of triangles) {
            let newVertices = [];
            const u = [tri[1][0] - tri[0][0], tri[1][1] - tri[0][1], tri[1][2] - tri[0][2]];
            const v = [tri[2][0] - tri[0][0], tri[2][1] - tri[0][1], tri[2][2] - tri[0][2]];

            for (let coef of baryCoefs) {
                newVertices.push(ssub(coef[1], ssub(coef[0], tri[0], u), v));
            }

            for (let trindex of trindices) {
                newTriangles.push([
                    newVertices[trindex[0]], newVertices[trindex[1]], newVertices[trindex[2]]
                ]);
            }
        }
        triangles = newTriangles.slice();
    }
    
    // normalize vertices
    for (let i=0; i<triangles.length; i++) {
        triangles[i][0] = snormal(triangles[i][0]);
        triangles[i][1] = snormal(triangles[i][1]);
        triangles[i][2] = snormal(triangles[i][2]);
    }
    return triangles;

}

function icosphere_flat(subdivisions=0) {
    const icosphere_tris = icosphere(subdivisions);
    const verts = [];
    for (const tri of icosphere_tris) {
        for (const vert of tri) {
            verts.push(vert);
        }
    }
    return verts;
}


module.exports = {
    snormal, ssub, sscale, // remove these later (Fix this hacky garbage)
    icosphere,
};