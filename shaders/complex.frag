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

vec2 asinC(vec2 z) {
    return divC(lnC(addC(multC(z, i), sqrtC(subC(vec2(1., 0.), multC(z, z))))), i);
}

vec2 acosC(vec2 z) {
    return divC(lnC(addC(z, sqrtC(subC(multC(z, z), vec2(1., 0.))))), i);
}

vec2 atanC(vec2 z) {
    return divC(lnC(divC(subC(i, z), addC(i, z))), vec2(0., 2.));
}

vec2 sinhC(vec2 z) {
    return scaleC( subC( expC(z), expC(scaleC(z, -1.)) ), 0.5 );
}

vec2 coshC(vec2 z) {
    return scaleC( addC( expC(z), expC(scaleC(z, -1.)) ), 0.5 );
}

vec2 tanhC(vec2 z) {
    return divC(sinhC(z), coshC(z));
}

vec2 atanhC(vec2 z) {
    return scaleC(lnC(divC(addC(vec2(1., 0.), z), subC(vec2(1., 0.), z))), 0.5);
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

vec2 GammaMainC(vec2 z) {
    z = subC(z, vec2(1., 0.));
    vec2 x = vec2(pValues[0], 0);
    for (int j=1; j<9; j++) {
        // there is quite literally no way to NOT hard code j<9 in the above loop
        // as in GLSL ES 2.0 you can't get the length of an array and you can't compare
        // the loop variable to a non-constant value (like a uniform)
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
    return addC(multC( subC(vec2(1., 0.), t), z ), multC(w, t) );
}

vec2 fracC(vec2 z) {
    return vec2(
        fract(z.x), fract(z.y)
    );
}

vec2 clampC(vec2 z, vec2 min_, vec2 max_) {
    float norm = normC(z).x;
    if (!(min_.x <= norm && norm <= max_.x)) {
        return scaleC(z, clamp(norm, min_.x, max_.x) / norm);
    }
    return z;
}

// SC stuff

vec2 nCrC(vec2 n, vec2 k) {
    return divC(GammaC(n + vec2(1., 0.)), multC( GammaC(k + vec2(1., 0.)), GammaC(n - k + vec2(1., 0.)) ));
}

vec2 scC(vec2 z, vec2 pC) {
    float p = pC.x;

    // compute scale factor (inverse radius of resulting polygon)
    float scaleFactor = p / betaC(vec2(1.0 / p, 0.), vec2(1.0 - 2.0 / p, 0.)).x;

    // compute integral
    float step = 1.0 / 100.;
    vec2 integral = vec2(0.0, 0.0);
    for (int i=0; i<101; i++) {
        float t = step * float(i);
        integral += step * multC( z, powC( vec2(1.0, 0.0) - powC( t * z, vec2(p, 0) ), vec2(-2.0 / p, 0) ) );
    }

    return scaleFactor * integral;
}

float scInvRadiusScale(float p) {
  float scaleFactor = 1.0 / cos(pi / p);
  return scaleFactor * ( 1.00508636015 - 0.673153865829 * pow(p, -1.55599992684) );
}

vec2 inverseSCC(vec2 z, vec2 pC) {
    float p = pC.x;

    // from what I can tell, there is no way to declare an array in GLSL ES 2.0,
    // which is what I did in the original version of this code (written in es 3.0) for the series coefficients
    // so this has to suffice
    float Cn1 = scaleC(nCrC(vec2(1. - 1.0 + 2.0 / p, 0.), vec2(1., 0.)), 1. / (1. + 1. + p)).x;
    float Cn2 = scaleC(nCrC(vec2(2. - 1.0 + 2.0 / p, 0.), vec2(2., 0.)), 1. / (1. + 2. + p)).x;
    float Cn3 = scaleC(nCrC(vec2(3. - 1.0 + 2.0 / p, 0.), vec2(3., 0.)), 1. / (1. + 3. + p)).x;
    float Cn4 = scaleC(nCrC(vec2(4. - 1.0 + 2.0 / p, 0.), vec2(4., 0.)), 1. / (1. + 4. + p)).x;
    float Cn5 = scaleC(nCrC(vec2(5. - 1.0 + 2.0 / p, 0.), vec2(5., 0.)), 1. / (1. + 5. + p)).x;

    float C = betaC(vec2(1.0 / p, 0.), vec2(1.0 - 2.0 / p, 0.)).x / p;

    // taylor coefficients
    float T1 = -Cn1;
    float T2 = -Cn2 + (p + 1.0) * pow(Cn1, 2.0);
    float T3 = -Cn3 + (3.0 * p + 2.0) * (Cn1 * Cn2 - (p + 1.0) / 2.0 * pow(Cn1, 3.0));
    float T4 = -Cn4 + (2.0 * p + 1.0) * (2.0 * Cn1 * Cn3 + pow(Cn2, 2.0) - (4.0 * p + 3.0) * (pow(Cn1, 2.0) * Cn2 - (p + 1.0) / 3.0 * pow(Cn1, 4.0)));
    float T5 = -Cn5 + (5.0 * p + 2.0) * (Cn1 * Cn4 + Cn2 * Cn3 + (5.0 * p + 3.0) * (-0.5 * pow(Cn1, 2.0) * Cn3 - 0.5 * Cn1 * pow(Cn2, 2.0) + (5.0 * p + 4.0) * ((pow(Cn1, 3.0) * Cn2) / 6.0 - ((p + 1.0) * pow(Cn1, 5.0)) / 24.0 )));

    z = scInvRadiusScale(p) * z;
    vec2 h = powC(z, vec2(p, 0.0));
    vec2 base = vec2(1.0, 0.0) + T1 * h + T2 * powC(h, vec2(2.0, 0.0)) + T3 * powC(h, vec2(3.0, 0.0)) + T4 * powC(h, vec2(4.0, 0.0)) + divC( T5 * powC(h, vec2(5.0, 0.0)), vec2(1.0, 0.0) + h / pow(C, p) );

    return multC(z, base);
}

vec2 planeToPC(vec2 z, vec2 pC) {
    vec2 squished = scaleC(z, 2. / (1. + exp(-normC(z).x)) - 1.);
    return scC(squished, pC);
}

vec2 pToPlaneC(vec2 z, vec2 pC) {
    z = inverseSCC(z, pC);
    float norm = normC(z).x;
    vec2 stretched = scaleC(z, log((1. + norm) / (1. - norm)));
    return stretched;
}

// end SC stuff

vec2 squeezeC(vec2 z, vec2 coverage, vec2 squeezeLength) {
    float norm = normC(z).x;
    float alpha = sqrt(atanhC(vec2(clamp(coverage.x, 0., 1.), 0.)).x / squeezeLength.x);
    float ratio = tanhC(vec2(alpha * alpha * norm, 0.)).x / norm;
    return scaleC(z, ratio);
}

vec2 floorC(vec2 z) {
    return vec2(floor(z.x), floor(z.y));
}

vec2 ceilC(vec2 z) {
    return vec2(ceil(z.x), floor(z.y));
}

vec2 modC(vec2 z, vec2 b) {
    return subC(z, multC(b, floorC(divC(z, b))));
}

bool fIsInvalid(float x) {
    return !(x <= 0. || 0. <= x) || (abs(x) > 1000000.0);
}

bool isInvalid(vec2 z) {
    return (fIsInvalid(z.x) || fIsInvalid(z.y));
}