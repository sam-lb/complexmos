


class rVector {

    constructor(x, y) {
        this.x = x;
        this.y = y;
    }

    dot(v) {
        return this.x * v.x + this.y * v.y;
    }

    magSq() {
        return this.x * this.x + this.y * this.y;
    }

    mag() {
        return Math.sqrt(this.magSq());
    }

    scale(k) {
        return new rVector(this.x * k, this.y * k);
    }

    proj(v) {
        return v.scale(this.dot(v) / v.magSq());
    }

    unit() {
        return this.scale(1 / this.mag());
    }

    angleBetween(v) {
        return Math.acos(this.dot(v).scale(this.mag() * v.mag()));
    }

    perp() {
        return new rVector(-this.y, this.x);
    }

    reflect(pointOnMirror, mirrorNormal) {
        return pointOnMirror.add(this.sub(pointOnMirror).proj(mirrorNormal)).scale(2).sub(this);
    }

    add(v) {
        return new rVector(this.x + v.x, this.y + v.y);
    }

    sub(v) {
        return new rVector(this.x - v.x, this.y - v.y);
    }

}


function rvec(x, y) {
    return new rVector(x, y);
}


module.exports = {
    rVector, rvec,
};