



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

function snorm(X) {
    return Math.sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
}

function sdist(X, Y) {
    return snorm(ssub(1, X, Y));
}

function sdot(X, Y) {
    return X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2];
}

function sangle(X, Y) {
    return Math.acos(sdot(X, Y) / (snorm(X) * snorm(Y)));
}

function scross(X, Y) {
    return [
        X[1] * Y[2] - X[2] * Y[1],
        X[2] * Y[0] - X[0] * Y[2],
        X[0] * Y[1] - X[1] * Y[0]
    ];
}

function sunormal(X, Y) {
    const norm = scross(X, Y);
    return snormal(norm);
}


// const REF_NORMAL = [0, 0, 1];
// function sortClockwise(points) {
//     const ref = ssub(1, points[1], points[0]);
    
//     let results = points.slice(1, points.length).sort((x, y) => {
//         return sangle(ref, ssub(1, x, points[0])) - sangle(ref, ssub(1, y, points[0]));
//     });
//     results.unshift(points[0]);

//     const normal = sunormal(ssub(1, results[1], results[0]), ssub(1, results[2], results[0]));
//     if (sangle(normal, REF_NORMAL) < Math.PI / 2) {
//         return results.reverse();
//     }
//     return results;
// }

// function sortClockwise(points) {
//     const normal = scross(ssub(1, points[1], points[0]), ssub(1, points[2], points[0]));
//     if (snorm(ssub(-1, points[0], normal)) < snorm(points[0])) {
//         return [points[0], points[2], points[1]];
//     }
//     return points;
// }

function scenter(points) {
    let xTotal = 0, yTotal = 0, zTotal = 0;
    for (let i=0; i<points.length; i++) {
        xTotal += points[i][0];
        yTotal += points[i][1];
        zTotal += points[i][2];
    }
    xTotal /= points.length;
    yTotal /= points.length;
    zTotal /= points.length;
    return [xTotal, yTotal, zTotal];
}

function sortClockwise(points) {
    const mid = scenter(points);
    const vec1 = ssub(1, points[0], mid);
    const vec2 = ssub(1, points[1], mid);
    const normal = scross(vec1, vec2);

    if (sdot(normal, mid) > 0) {
        return points;
    } else {
        return [points[0], points[2], points[1]];
    }
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

    for (let i=0; i<triangles.length; i++) {
        triangles[i] = sortClockwise(triangles[i]);
    }

    return triangles;

}

function icosphere_flat(subdivisions=0) {
    return flatten(icosphere(subdivisions));
}

function icosphere_flat_lopsided(subdivisions_top=0, subdivisions_bottom=0) {
    const top = icosphere(subdivisions_top).filter(face => {
        return scenter(face)[2] < 0;
    });
    const bottom = icosphere(subdivisions_bottom).filter(face => {
        return scenter(face)[2] >= 0;
    });
    return flatten(Array.prototype.concat(top, bottom));
}

function flatten(triangles) {
    const verts = [];
    for (const tri of triangles) {
        for (const vert of tri) {
            verts.push(vert);
        }
    }
    return verts;
}


module.exports = {
    snormal, ssub, sscale, // remove these later (Fix this hacky garbage)
    icosphere, icosphere_flat, icosphere_flat_lopsided
};