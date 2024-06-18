


class TokenType {

    static LEFT_PAREN = 1;
    static RIGHT_PAREN = 2;
    static COMMA = 3;
    static ASSIGN = 4;
    static PLUS = 5;
    static MINUS = 6;
    static ASTERISK = 7;
    static SLASH = 8;
    static CARET = 9;
    static COLON = 10;
    static NAME = 11;
    static NUMBER = 12;
    static EOF = 13;

    constructor(value) {
        this.value = value;
    }

    static toStr(tokenType) {
        const values = [
            "(", ")", ",", "=", "+", "-", "*", "/", "^", ":", "", "", "",
        ];
        return values[tokenType - 1];
    }

}

module.exports = {
    TokenType,
};