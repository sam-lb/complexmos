vec3 hsvToRgb(vec3 c){
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 getColorAtIndex(int index) {
    // I don't want to pass the color values as a texture. So this is the best solution, because you can't index arrays with non-constants.
    if (index == 0) {
        return vec3(gradR[0], gradG[0], gradB[0]);
    } else if (index == 1) {
        return vec3(gradR[1], gradG[1], gradB[1]);
    } else if (index == 2) {
        return vec3(gradR[2], gradG[2], gradB[2]);
    } else if (index == 3) {
        return vec3(gradR[3], gradG[3], gradB[3]);
    } else if (index == 4) {
        return vec3(gradR[4], gradG[4], gradB[4]);
    } else if (index == 5) {
        return vec3(gradR[5], gradG[5], gradB[5]);
    } else {
        return vec3(0., 0., 0.);
    }
}

vec3 getColorGradient(vec2 outp) {
    float angle = atan(outp.y, outp.x);
    float nm = normC(outp).x;
    float trm = .25 + .75 * floor((2. / pi * atan(sqrt(nm))) / 0.2) * 0.2;
    float sector = tpi / 6.;
    float i0 = mod(floor(angle / sector), 6.) + 1.;
    float i1 = mod(i0, 6.) + 1.;
    vec3 firstColor = getColorAtIndex(int(i0) - 1);
    vec3 secondColor = getColorAtIndex(int(i1) - 1);
    float t0 = (angle - floor(angle / sector)) / sector;

    return trm * ((1. - t0) * firstColor + t0 * secondColor);
}

vec3 getColorDefault(vec2 outp) {
    float nm = normC(outp).x;
    float trm = .25 + .75 * floor((2. / pi * atan(sqrt(nm))) / 0.2) * 0.2;
    vec3 col = vec3(mod(atan(outp.y, outp.x) + tpi,  tpi) / tpi, 1., trm);
    return hsvToRgb(col);
}

vec3 getColorTextureRepeat(vec2 outp) {
    vec2 texCoord = vec2(
        outp.x - floor(outp.x),
        1. - (outp.y - floor(outp.y))
    );
    return texture2D(texture, texCoord).xyz;
}

vec3 getColorTextureStretch(vec2 outp) {
    float squareNorm = max(abs(outp.x), abs(outp.y));
    vec2 texCoord = 1. / pi * atan(squareNorm) * outp / squareNorm + vec2(.5, .5);
    texCoord = vec2(texCoord.x, 1. - texCoord.y);
    return texture2D(texture, texCoord).xyz;
}