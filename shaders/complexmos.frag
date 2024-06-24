precision highp float;
uniform float width, height;

/* complex lib */

const float pi = 3.1415926535897;
const float tpi = 2.0 * pi;
const vec2 i = vec2(0., 1.);
const float EPSILON = 0.0000001;

const float[9] pValues = float[](
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
);

vec2 conjC(vec2 z) {
    return vec2(z.x, -z.y);
}

vec2 normC(vec2 z) {
    return vec2(sqrt(z.x * z.x + z.y * z.y), 0.0);
}

vec2 normSqC(vec2 z) {
    return vec2(z.x * z.x + z.y * z.y, 0.0);
}

vec2 addC(vec2 z, vec2 w) {
    // obviously such functions are not stricly necessary but it's for ease of conversion from AST to glsl
    // maybe it's optimized away by the compiler anyway
    return vec2(z.x + w.x, z.y + w.y);
}

vec2 subC(vec2 z, vec2 w) {
    return vec2(z.x - w.x, z.y - w.y);
}

vec2 multC(vec2 z, vec2 w) {
    return vec2(z.x * w.x - z.y * w.y, z.x * w.y + z.y * w.x);
}

vec2 scaleC(vec2 z, float k) {
    return z * k;
}

vec2 invC(vec2 z) {
    return scaleC(conjC(z), 1. / normSqC(z).x);
}

vec2 divC(vec2 z, vec2 w) {
    return multC(z, invC(w));
}

vec2 argC(vec2 z) {
    return vec2(mod(atan(z.y, z.x) + tpi, tpi), 0.0);
}

vec2 expC(vec2 z) {
    float mag = exp(z.x);
    return vec2(cos(z.y), sin(z.y)) * mag;
}

vec2 lnC(vec2 z) {
    return vec2(log(normC(z).x), argC(z).x);
}

vec2 sqrtC(vec2 z) {
    float normSqrt = sqrt(normC(z).x);
    float halfArg = 0.5 * argC(z).x;
    return vec2(cos(halfArg), sin(halfArg)) * normSqrt;
}

vec2 sinC(vec2 z) {
    vec2 rotated = multC(z, i);
    return divC(subC(expC(rotated), expC(scaleC(rotated, -1.))), scaleC(i, 2.));
}

vec2 cosC(vec2 z) {
    vec2 rotated = multC(z, i);
    return scaleC(addC(expC(rotated), expC(scaleC(rotated, -1.))), 0.5);
}

vec2 tanC(vec2 z) {
    return divC(sinC(z), cosC(z));
}

vec2 sinhC(vec2 z) {
    return scaleC( subC( expC(z), expC(scaleC(z, -1.)) ) );
}

vec2 coshC(vec2 z) {
    return scaleC( addC( expC(z), expC(scaleC(z, -1.)) ) );
}

vec2 tanhC(vec2 z) {
    return divC(sinhC(z), coshC(z));
}

vec2 reC(vec2 z) {
    return vec2(z.x, 0.0);
}

vec2 imC(vec2 z) {
    return vec2(z.y, 0.0);
}

vec2 powC(vec2 z, vec2 w) {
    float nSq = normSqC(z).x;

    if (nSq < EPSILON) {
        if (normSqC(subC(w, vec2(1., 0.))).x < EPSILON) {
            return vec2(1., 0.);
        } else {
            return vec2(0., 0.);
        }
    }

    float subAng = atan(z.y, z.x);
    float ang = 0.5 * w.y * log(nSq) + w.x * subAng;
    float nm = exp(-w.y * subAng) * pow(nSq, 0.5 * w.x);
    return vec2(cos(ang), sin(ang)) * nm;
}

z = Complex.sub(z, complex(1, 0)); // account for stupid shift by 1
let x = complex(pValues[0], 0);
for (let i=1; i<pValues.length; i++) {
    x = Complex.add(x, Complex.div( complex(pValues[i], 0.0), Complex.add(z, complex(i, 0)) ));
}
const t = Complex.add(z, complex(7.5, 0)); // g=7, g+0.5
return Complex.pow(t, z.add(complex(0.5, 0))).mult(t.scale(-1).exp()).mult(x).scale(Math.sqrt(2 * Math.PI));

vec2 GammaMainC(vec2 z) {
    z = subC(z, vec2(1., 0.));
    vec2 x = vec2(pValues[0], 0);
    for (int j=1; j<pValues.length(); j++) {
        x = addC(x, divC( vec2(pValues[j], 0.0), addC(z, vec2(float(j), 0.0)) ));
    }
    vec2 t = addC(z, vec2(7.5, 0.0)); // g=7, g+0.5
    return scaleC(multC(multC(powC(t, addC(z, vec2(0.5, 0))), expC(scaleC(t, -1.))), x), sqrt(tpi));
}

vec2 GammaC(vec2 z) {
    if (z.x < 0.5) {
        // gamma(1-z)gamma(z) = pi / sin(pi * z)
        return divC(vec2(pi, 0.0), multC( sinC(scaleC(z, pi)), GammaMainC(vec2(1. - z.x, -z.y)) ));
    } else {
        return GammaMainC(z);
    }
}

vec2 betaC(vec2 z, vec2 w) {
    return divC( multC(GammaC(z), GammaC(w)), GammaC(addC(z, w)) );
}

vec2 minC(vec2 z, vec2 w) {
    if (normSqC(z).x < normSqC(w).x) {
        return z;
    }
    return w;
}

vec2 maxC(vec2 z, vec2 w) {
    if (normSqC(z).x < normSqC(w).x) {
        return w;
    }
    return z;
}

vec2 lerpC(vec2 z, vec2 w, vec2 t) {
    addC(multC( subC(vec2(1., 0.), t), z ), multC(w, t) );
}

/* end complex lib */

void main() {
    gl_FragColor = vec4(gl_FragCoord.x/width, gl_FragCoord.y/height, 0.0, 1.0);
}