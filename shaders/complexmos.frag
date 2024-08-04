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

//IMPORT_COLORING

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