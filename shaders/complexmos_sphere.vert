attribute vec3 position;

uniform float width, height;
uniform vec3 row1;
uniform vec3 row2;
uniform vec3 row3;

varying vec3 outPos;

vec3 transform(vec3 P) {
    return vec3(
        row1.x * P.x + row1.y * P.y + row1.z * P.z,
        row2.x * P.x + row2.y * P.y + row2.z * P.z,
        row3.x * P.x + row3.y * P.y + row3.z * P.z
    );
}

void main() {
    vec3 transformed = transform(position);
    outPos = position;
    gl_Position = vec4(transformed.x / (-1.5 + transformed.y) * height / width * 0.9, transformed.z / (-1.5 + transformed.y) * 0.9, -transformed.y, 1.0);
}