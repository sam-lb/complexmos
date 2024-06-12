


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

}


function rvec(x, y) {
    return new rVector(x, y);
}