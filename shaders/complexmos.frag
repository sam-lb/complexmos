precision highp float;
uniform float width, height;
uniform sampler2D texture;

const float pi = 3.1415926535897;
const float tpi = 2.0 * pi;
const vec2 e = vec2(2.71828182846, 0.);
const vec2 piC = vec2(pi, 0.);
const vec2 tpiC = vec2(tpi, 0.);
const vec2 i = vec2(0., 1.);
const float EPSILON = 0.0000001;

uniform float pValues[9]; // GLSL ES 2.0 sucks. There's no way to initialize an array
uniform float gradR[6]; // GLSL ES 2.0 sucks. can't index into an array with non-constant. Regl can't pass in vec3 arrays as uniforms
uniform float gradG[6];
uniform float gradB[6];
uniform vec2 xBounds;
uniform vec2 yBounds;

//IMPORT_COMPLEX

/* end complex lib */

//REPLACE_BEGIN
vec2 udf_f(vec2 z) {
    return vec2(1., 0.);
}
//REPLACE_END

float round(float x) {
    if (fract(x) > 0.5) {
        return ceil(x);
    }
    return floor(x);
}

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

vec3 monokai(float angle) {
    float sector = tpi / 6.;
    float i0 = mod(floor(angle / sector), 6.) + 1.;
    float i1 = mod(i0, 6.) + 1.;
    vec3 firstColor = getColorAtIndex(int(i0) - 1);
    vec3 secondColor = getColorAtIndex(int(i1) - 1);
    // vec4 firstColor = texture2D(gradColors, vec2((i0 - 1. + 0.5) / 6., 0.));
    // vec4 secondColor = texture2D(gradColors, vec2((i1 - 1. + 0.5) / 6.), 0.);
    // float firstColor = texture2D(gradColors, vec2(0., 0.));
    // float secondColor = texture2D(gradColors, vec2(0., 0.));
    float t0 = (angle - floor(angle / sector)) / sector;

    return (1. - t0) * firstColor + t0 * secondColor;
}

void main() {

    float xUnits = xBounds.y - xBounds.x;
    float yUnits = yBounds.y - yBounds.x;
    float xCenter = xBounds.x + 0.5 * xUnits;
    float yCenter = yBounds.x + 0.5 * yUnits;

    vec2 uv = vec2(gl_FragCoord.x / width, gl_FragCoord.y / height);
    float aspect = yUnits / xUnits;
    vec2 z = ((uv - vec2(0.5, 0.5)) * xUnits + vec2(xCenter, yCenter)) * vec2(1, aspect);

//DISPLAY_REPLACE_BEGIN
    vec2 outp = udf_f(z);
//DISPLAY_REPLACE_END

    vec3 col;

    // float squareNorm = max(abs(outp.x), abs(outp.y));
    // vec2 texCoord = 1. / pi * atan(squareNorm) * outp / squareNorm + vec2(.5, .5);
    // texCoord = vec2(texCoord.x, 1. - texCoord.y);

    // vec2 texCoord = vec2(
    //     outp.x - floor(outp.x),
    //     1. - (outp.y - floor(outp.y))
    // );

    // vec4 col = texture2D(texture, texCoord);

    float nm = normC(outp).x;
    float trm = .25 + .75 * floor((2. / pi * atan(sqrt(nm))) / 0.2) * 0.2;
    // col = vec3(mod(atan(outp.y, outp.x) + tpi,  tpi) / tpi, 1., trm);
    // col = hsvToRgb(col);
    col = monokai(atan(outp.y, outp.x));

    float tolerance = 0.01;
    if (abs(z.x - round(z.x)) < tolerance || abs(z.y - round(z.y)) < tolerance) {
        col = vec3((0.25 + col.x) / 2., (0.25 + col.y) / 2., (0.25 + col.z) / 2.);
        // col = vec4((0.25 + col.x) / 2., (0.25 + col.y) / 2., (0.25 + col.z) / 2., col.w);
    }

    // gl_FragColor = col;
    gl_FragColor = vec4(col, 1.0);
}