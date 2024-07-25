precision highp float;
uniform float width, height;

varying vec2 outPos;

uniform float pValues[9]; // GLSL ES 2.0 sucks. There's no way to initialize an array
uniform vec2 xBounds;
uniform vec2 yBounds;

uniform float alpha;
uniform float beta;

const float pi = 3.1415926535897;
const float tpi = 2.0 * pi;
const vec2 e = vec2(2.71828182846, 0.);
const vec2 piC = vec2(pi, 0.);
const vec2 tpiC = vec2(tpi, 0.);
const vec2 i = vec2(0., 1.);
const float EPSILON = 0.0000001;

//IMPORT_COMPLEX

//REPLACE_BEGIN
vec2 udf_f(vec2 z) {
    return vec2(1., 0.);
}
//REPLACE_END

vec2 stereoProject(vec3 P) {
    float denom = 1. + P.z;
    return vec2(P.x / denom, P.y / denom);
}

vec3 hsvToRgb(vec3 c){
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec2 mandelbrot(vec2 z) {
    vec2 z0 = vec2(0., 0.);
    for (int i=0; i<50; i++) {
        z0 = addC(multC(z0, z0), z);
    }
    return z0;
}

float round(float x) {
    if (fract(x) > 0.5) {
        return ceil(x);
    }
    return floor(x);
}

void main() {
    float xUnits = xBounds.y - xBounds.x;
    float yUnits = yBounds.y - yBounds.x;
    float xCenter = xBounds.x + 0.5 * xUnits;
    float yCenter = yBounds.x + 0.5 * yUnits;

    vec2 uv = vec2(outPos.x, outPos.y);
    float aspect = yUnits / xUnits;
    vec2 z = (uv * xUnits + vec2(xCenter, yCenter)) * vec2(1, aspect);

    vec3 col;
    
//DISPLAY_REPLACE_BEGIN
    vec2 outp = udf_f(z);
//DISPLAY_REPLACE_END

    if (isInvalid(outp)) {
        col = vec3(1., 1., 1.);
    } else {
        float nm = normC(outp).x;
        float trm = .25 + .75 * floor((2. / pi * atan(sqrt(nm))) / 0.2) * 0.2;
        col = vec3(mod(atan(outp.y, outp.x) + tpi,  tpi) / tpi, 1., trm);
        col = hsvToRgb(col);
    }

    float tolerance = 0.01;
    if (abs(z.x - round(z.x)) < tolerance || abs(z.y - round(z.y)) < tolerance) {
        col = vec3((0.25 + col.x) / 2., (0.25 + col.y) / 2., (0.25 + col.z) / 2.);
    }

    gl_FragColor = vec4(col, 1.0);
}