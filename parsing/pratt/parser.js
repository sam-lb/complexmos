


const { tracker } = require("../errors.js");
const { TokenType } = require("./tokentype.js");
const { Precedence } = require("./precedence.js");
const { Token } = require("./token.js");


class Parser {

    constructor(tokens) {
        this.mTokens = tokens;
        this.mPrefixParselets = {};
        this.mInfixParselets = {};
        this.index = 0;
    }

    registerPrefix(tokenType, parselet) {
        this.mPrefixParselets[tokenType] = parselet;
    }

    registerInfix(tokenType, parselet) {
        this.mInfixParselets[tokenType] = parselet;
    }

    parseExpression(precedence=Precedence.LOWEST) {
        let token = this.consume();
        if (tracker.hasError) return;
        if (!Object.keys(this.mPrefixParselets).includes(token.mtype.toString())) {
            if (token.mtype === TokenType.EOF) {
                tracker.error(`Unexpected EOF while parsing (1)`);
            } else {
                tracker.error(`Invalid syntax`);
            }
            return;
        }

        let left = this.mPrefixParselets[token.mtype].parse(this, token);
        if (tracker.hasError) return;

        while (precedence < this.getPrecedence()) {
            token = this.consume();
            if (Object.keys(this.mInfixParselets).includes(token.mtype.toString())) {
                left = this.mInfixParselets[token.mtype].parse(this, left, token);
                if (tracker.hasError) return;
            }
        }

        return left;
    }

    consume(expected=null) {
        const token = this.lookAhead(0);
        if (expected !== null && token.mtype === TokenType.EOF) {
            tracker.error(`Unexpected EOF while parsing`);
            return;
        } else if (expected !== null && token.mtype !== expected) {
            tracker.error(`Unexpected EOF while parsing`);
            return;
        }
        this.index++;
        return token;
    }

    lookAhead(distance) {
        if (this.index + distance < this.mTokens.length) {
            return this.mTokens[this.index + distance];
        } else {
            return new Token(TokenType.EOF, null);
        }
    }

    peek(expected) {
        const token = this.lookAhead(0);
        return (token.mtype === expected);
    }

    getPrecedence() {
        const currentToken = this.lookAhead(0);
        if (Object.keys(this.mInfixParselets).includes(currentToken.mtype.toString())) {
            return this.mInfixParselets[currentToken.mtype].getPrecedence();
        }
        return Precedence.LOWEST;
    }

}


module.exports = {
    Parser,
};