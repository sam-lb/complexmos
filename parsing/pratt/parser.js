


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
        if (!Object.keys(this.mPrefixParselets).includes(token.mtype.toString())) {
            console.error(`couldn't parse token ${token.toString()}`);
        }

        let left = this.mPrefixParselets[token.mtype].parse(this, token);

        while (precedence < this.getPrecedence()) {
            token = this.consume();
            if (Object.keys(this.mInfixParselets).includes(token.mtype.toString())) {
                left = this.mInfixParselets[token.mtype].parse(this, left, token);
            }
        }

        return left;
    }

    consume(expected=null) {
        const token = this.mTokens[this.index];
        if (expected !== null && token.mtype !== expected) {
            console.error(`consumed token ${token.mtype} does not match expected ${expected}`);
        }
        this.index++;
        return token;
    }

    lookAhead(distance) {
        if (this.index + distance < this.mTokens.length) {
            return this.mTokens[this.index + distance];
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

