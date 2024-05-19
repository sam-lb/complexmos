

function stereographic(z) {
    /**
     * computes the stereographic projection of the 3d point z onto the xy plane
     */
    return complex(
        z[0] / (1 - z[2]),
        z[1] / (1 - z[2]),
    );
}

function inverseStereoProject(z) {
    /**
     * computes the inverse stereographic projection of a 2d point onto the unit spheree
     */
    const normSq = z.normSq();
        return matrix([
            (2 * z.re) / (1 + normSq),
            (2 * z.im) / (1 + normSq),
            (normSq - 1) / (1 + normSq),
        ]).transpose();
}

function perspectiveProject(Z, orbit, fov) {
    /**
     * computes the perspective projection of the 3d point Z onto the xz plane
     * given the orbit radius and fov values
     */
    const values = Z.getColumn(0);
    const x = values[0], y = values[1], z = values[2];
    const denom = (orbit + fov * y);
    return complex(
        x / denom, z / denom,
    );
}