attribute vec2 position;
uniform float angle, width, height;
void main() {
    float aspect = width / height;
    gl_Position = vec4(
        (cos(angle) * position.x - sin(angle) * position.y),
        aspect * (sin(angle) * position.x + cos(angle) * position.y),
        0.0,
        1.0
    );
}