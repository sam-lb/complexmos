precision highp float;
uniform float width, height;

void main() {
    gl_FragColor = vec4(gl_FragCoord.x/width, gl_FragCoord.y/height, 0.0, 1.0);
}