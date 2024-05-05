


function parameterizePoints(points) {
    /*
    Return a parameteric function consisting of straight lines connecting points
    Note: does NOT necessarily produce a constant time parameterization
    */
    const N = points.length;
    const segmentLength = 1 / (N - 1);
    return (t) => {
        const index = (t === 1) ? N - 2 : Math.floor(t / segmentLength);
        return Euclid.lerp(points[index], points[index+1], (t - (t % segmentLength)) / segmentLength);
    };
}



function fourierCoefficients(f, N) {
    /*
    Calculate the fourier coefficients for the complex function f 
    parameterized on [0, 1]. 
    Coefficients from -N, -N+1, ... , N-1, N
    */
    const coefs = [];
    for (let n=-N; n<=N; n++) {
        const F = (t) => {
            return Complex.mult(f(t), Complex.exp( complex(0, -2 * Math.PI * n) ));
        };
        coefs.push(
            integrateOverParameter(F, 0, 1)
        );
    }
    return coefs;
}


function fourierSeries(f, N) {
    const coefs = fourierCoefficients(f, N);
    return (t) => {
        const result = complex(0, 0);
        for (let n=-N; n<=N; n++) {
            const coef = coefs[n + N];
            result.iadd(
                Complex.mult(coef, Complex.exp( complex(0, 2 * Math.PI * n) ))
            );
        }
        return result;
    };
}