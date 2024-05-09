



function icosphere(subdivisions=0) {

    const scale = 1 / (2 * Math.sin(2 * Math.PI / 5));
    const phi = scale * (1 + Math.sqrt(5)) / 2;

    /**
     * [-s,s,-s,s,0,0,0,0,phi,phi,-phi,-phi]
     * [phi,phi,-phi,-phi,-s,s,-s,s,0,0,0,0]
     * [0,0,0,0,phi,phi,-phi,-phi,-s,s,-s,s]
     */

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
    
    return triangles;

}