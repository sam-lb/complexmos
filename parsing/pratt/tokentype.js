


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
    static NAME = 10;
    static NUMBER = 11;
    static EOF = 12;

    constructor(value) {
        this.value = value;
    }

    static toStr(tokenType) {
        const values = [
            "(", ")", ",", "=", "+", "-", "*", "/", "^", "", "", "",
        ];
        return values[tokenType - 1];
    }

}