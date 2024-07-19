attribute vec3 position;

uniform float width, height;
uniform vec3 row1;
uniform vec3 row2;
uniform vec3 row3;

uniform float pValues[9]; // GLSL ES 2.0 sucks. There's no way to initialize an array
uniform vec2 xBounds;
uniform vec2 yBounds;

const float pi = 3.1415926535897;
const float tpi = 2.0 * pi;
const vec2 e = vec2(2.71828182846, 0.);
const vec2 piC = vec2(pi, 0.);
const vec2 tpiC = vec2(tpi, 0.);
const vec2 i = vec2(0., 1.);
const float EPSILON = 0.0000001;

varying vec2 outPos;

//IMPORT_COMPLEX

//REPLACE_BEGIN
vec2 udf_f(vec2 z) {
    return vec2(0., 0.);
}
//REPLACE_END

vec3 transform(vec3 P) {
    return vec3(
        row1.x * P.x + row1.y * P.y + row1.z * P.z,
        row2.x * P.x + row2.y * P.y + row2.z * P.z,
        row3.x * P.x + row3.y * P.y + row3.z * P.z
    );
}

void main() {
    outPos = vec2(position.x * height / width, position.y);
    
    float xUnits = xBounds.y - xBounds.x;
    float yUnits = yBounds.y - yBounds.x;
    float xCenter = xBounds.x + 0.5 * xUnits;
    float yCenter = yBounds.x + 0.5 * yUnits;

    vec2 uv = vec2(outPos.x, outPos.y);
    float aspect = yUnits / xUnits;
    vec2 z = (uv * xUnits + vec2(xCenter, yCenter)) * vec2(1, aspect);

    vec3 newPosition = vec3(position.x, position.y, -normC(udf_f(z)).x);

    vec3 transformed = transform(newPosition);   
    vec2 projected = vec2(transformed.x / (-1.5 + transformed.y) * 0.9, transformed.z / (-1.5 + transformed.y) * 0.9);
    gl_Position = vec4(projected, -transformed.y, 1.0);
}