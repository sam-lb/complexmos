

class Token {

    constructor(tokenType, text) {
        this.mtype = tokenType;
        this.text = text;
    }

    toString() {
        return `<${this.mtype}, ${this.text}>`;
    }

}


module.exports = {
    Token,
};