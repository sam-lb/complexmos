

function stereographic(z) {
    // z = [
    //     (z[0] + 1) / 2,
    //     (z[1] + 1) / 2,
    //     (z[2] + 1) / 2,
    // ];

    return complex(
        z[0] / (1 - z[2]),
        z[1] / (1 - z[2]),
    );
}