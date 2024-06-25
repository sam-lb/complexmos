attribute vec2 position;

void main() {
    gl_Position = vec4(0.99 * position.x, 0.99 * position.y, 0.0, 1.0);
}