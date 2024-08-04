precision highp float;
uniform float width, height;
uniform sampler2D texture;

varying vec3 outPos;

uniform float pValues[9]; // GLSL ES 2.0 sucks. There's no way to initialize an array
uniform vec2 xBounds;
uniform vec2 yBounds;

uniform float alpha;
uniform float beta;

uniform float gradR[6]; // GLSL ES 2.0 sucks. can't index into an array with non-constant. Regl can't pass in vec3 arrays as uniforms
uniform float gradG[6];
uniform float gradB[6];

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

//IMPORT_COLORING

vec2 stereoProject(vec3 P) {
    float denom = 1. + P.z;
    return vec2(P.x / denom, P.y / denom);
}

float round(float x) {
    if (fract(x) > 0.5) {
        return ceil(x);
    }
    return floor(x);
}

void main() {
    float xCenter = 0.5 * (xBounds.x + xBounds.y);
    float yCenter = 0.5 * (yBounds.x + yBounds.y);

    vec2 projected = stereoProject(vec3(outPos.x, outPos.y, outPos.z));
    float aspect = height / width;
    vec2 z = projected * aspect + vec2(xCenter, yCenter);
    
//DISPLAY_REPLACE_BEGIN
    vec2 outp = udf_f(z);
    vec3 col;
    if (isInvalid(outp)) {
        col = vec3(0., 0., 0.);
    } else {
        col = getColorDefault(outp);
    }
//DISPLAY_REPLACE_END

    float tolerance = 0.01;
    if (abs(z.x - round(z.x)) < tolerance || abs(z.y - round(z.y)) < tolerance) {
        col = vec3((0.25 + col.x) / 2., (0.25 + col.y) / 2., (0.25 + col.z) / 2.);
    }

    gl_FragColor = vec4(col, 1.0);
}