

const alnum = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";



class Lexer {

    constructor(mText) {
        this.mPunctuators = [];
        this.mText = mText;
        this.index = 0;
    }

    getTokenName() {
        let identifier = "";
        while (this.index < this.mText.length && alnum.includes(this.mText[this.index])) {
            identifier += this.mText[this.index];
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.NAME, identifier));
    }

    tokenize() {
        while (this.index < this.mText.length) {
            const char = this.mText[this.index];
            
            if (char === "(") {
                this.mPunctuators.push(new Token(TokenType.LEFT_PAREN, char));
            } else if (char === ")") {
                this.mPunctuators.push(new Token(TokenType.RIGHT_PAREN, char));
            } else if (char === ",") {
                this.mPunctuators.push(new Token(TokenType.COMMA, char));
            } else if (char === "=") {
                this.mPunctuators.push(new Token(TokenType.ASSIGN, char));
            } else if (char === "+") {
                this.mPunctuators.push(new Token(TokenType.PLUS, char));
            } else if (char === "-") {
                this.mPunctuators.push(new Token(TokenType.MINUS, char));
            } else if (char === "*") {
                this.mPunctuators.push(new Token(TokenType.ASTERISK, char));
            } else if (char === "/") {
                this.mPunctuators.push(new Token(TokenType.SLASH, char));
            } else if (char === "^") {
                this.mPunctuators.push(new Token(TokenType.CARET, char));
            } else if (alnum.includes(char)) {
                this.getTokenName();
                continue;
            } else {
                console.error(`bruh what is this character even doing in your input ${char}`);
                this.index++;
                continue;
            }
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.EOF, null));
    }

    getTokens() {
        return this.mPunctuators;
    }

}