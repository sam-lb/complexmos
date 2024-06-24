precision highp float;
uniform float width, height;

/* complex lib */

vec2 conj(vec2 z) {
    return vec2(z.x, -z.y);
}

/* end complex lib */

void main() {
    gl_FragColor = vec4(gl_FragCoord.x/width, gl_FragCoord.y/height, 0.0, 1.0);
}