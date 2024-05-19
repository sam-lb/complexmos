

function stereographic(z) {
    return complex(
        z[0] / (1 - z[2]),
        z[1] / (1 - z[2]),
    );
}

function inverseStereoProject(z) {
    const normSq = z.normSq();
        return matrix([
            (2 * z.re) / (1 + normSq),
            (2 * z.im) / (1 + normSq),
            (normSq - 1) / (1 + normSq),
        ]).transpose();
}

function perspectiveProject(Z, orbit, fov) {
    const values = Z.getColumn(0);
    const x = values[0], y = values[1], z = values[2];
    const denom = (orbit + fov * y);
    return complex(
        x / denom, z / denom,
    );
}