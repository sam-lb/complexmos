

const num = "0123456789";
const alnum = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";



class Lexer {

    constructor(mText, allowUnboundIdentifiers=false) {
        this.mPunctuators = [];
        this.mText = mText;
        this.index = 0;
        this.allowUnboundIdentifiers = allowUnboundIdentifiers;
        this.setScope({});
    }

    setScope(scope) {
        this.scope = scope;
        this.builtinLookup = new Trie(Object.keys(scope.builtin === undefined ? {} : scope.builtin));
        this.userGlobalLookup = new Trie(Object.keys(scope.userGlobal === undefined ? {} : scope.userGlobal));
    }

    setAllowUnboundIdentifiers(allowUnboundIdentifiers) {
        this.allowUnboundIdentifiers = allowUnboundIdentifiers;
    }

    setText(mText) {
        this.mPunctuators = [];
        this.mText = mText;
        this.index = 0;
    }

    getTokenName() {
        // consume characters until non-identifier character is encountered, add them to buffer
        let buffer = "";
        while (this.index < this.mText.length && alnum.includes(this.mText[this.index])) {
            buffer += this.mText[this.index];
            this.index++;
        }

        // split buffer into identifiers, add implicit multiplication where necessary
        const identifiers = [];
        while (buffer.length > 0) {
            let matchFound = false;
            let possibleIdentifier, i;
            for (i=buffer.length-1; i>=0; i--) {
                possibleIdentifier = buffer.slice(0, i+1);
                if (this.builtinLookup.containsKey(possibleIdentifier) || this.userGlobalLookup.containsKey(possibleIdentifier)) {
                    matchFound = true;
                    break;
                }
            }
            if (matchFound) {
                identifiers.push(possibleIdentifier);
                buffer = buffer.slice(i+1, buffer.length);
            } else {
                i++; // no break statement
                if (alnum.includes(possibleIdentifier[0]) && !num.includes(possibleIdentifier[0])) {
                    if (this.allowUnboundIdentifiers) {
                        identifiers.push(possibleIdentifier);
                        buffer = buffer.slice(i+1, buffer.length);
                    } else {
                        tracker.error(`Undefined identifier ${buffer}`);
                        this.kill();
                        return;
                    }
                } else {
                    tracker.error("unexpected number following identifier");
                    this.kill();
                    return;
                }
            }
        }
        for (let identifier of identifiers) {
            this.mPunctuators.push(new Token(TokenType.NAME, identifier));
            this.mPunctuators.push(new Token(TokenType.ASTERISK, "*"));
        }
        this.mPunctuators.pop();
    }

    getTokenNumber() {
        let number = "";
        let decimalEncounted = false;
        while (this.index < this.mText.length && (num.includes(this.mText[this.index]) || this.mText[this.index] === ".")) {
            if (this.mText[this.index] === ".") {
                if (decimalEncounted) {
                    tracker.error("two decimal points encountered in number");
                    this.kill();
                    return;
                }
                decimalEncounted = true;
            }
            number += this.mText[this.index];
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.NUMBER, number));
    }

    kill() {
        /** make tokenization stop */
        this.index = this.mText.length;
    }

    tokenize() {
        this.mPunctuators = [];
        this.index = 0;

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
            } else if (num.includes(char)) {
                this.getTokenNumber();
                continue;
            } else if (alnum.includes(char) || char === ".") {
                this.getTokenName();
                continue;
            } else {
                tracker.error(`bruh what is this character even doing in your input ${char}`);
                this.kill();
                return;
            }
            this.index++;
        }
        this.mPunctuators.push(new Token(TokenType.EOF, null));

        /**
         * Implicit multiplication loop
         * cases handled here (let x, y = identifiers and n = number ):
         * xn -> syntax error
         * nx -> n * x
         * n( -> n * (
         * )n -> ) * n
         * )x -> ) * x
         * x( -> x(         (assuming scope[x].isFunction === false)
         * case not handled here:
         * xy -> x * y. this case is handled in getTokenName()
         */
        const tokens = (this.mPunctuators.length > 0) ? [this.mPunctuators[0]] : [];
        for (let i=0; i<this.mPunctuators.length-1; i++) {
            const token = this.mPunctuators[i];
            const nextToken = this.mPunctuators[i+1];
            let needsMultiplication = false;

            if (token.mtype === TokenType.NAME && nextToken.mtype === TokenType.NUMBER) {
                // this case will theoretically never happen (it'll raise an error in getTokenName()), but it's here for completeness
                tracker.error("unexpected number following identifier");
                this.kill();
                return;
            } else if ( // obviously, you can reduce the number of comparisons in this condition, but it's more clear this way
                token.mtype === TokenType.NUMBER && nextToken.mtype === TokenType.NAME ||           // nx
                token.mtype === TokenType.NUMBER && nextToken.mtype === TokenType.LEFT_PAREN ||     // n(
                token.mtype === TokenType.RIGHT_PAREN && nextToken.mtype === TokenType.NUMBER ||    // )n
                token.mtype === TokenType.RIGHT_PAREN && nextToken.mtype === TokenType.NAME         // )x
            ) {
                needsMultiplication = true;
            } else if (token.mtype === TokenType.NAME && nextToken.mtype === TokenType.LEFT_PAREN) {
                if (this.builtinLookup.containsKey(token.text)) {
                    needsMultiplication = !this.scope.builtin[token.text].isFunction;
                } else if (this.userGlobalLookup.containsKey(token.text)) {
                    needsMultiplication = !this.scope.userGlobal[token.text].isFunction;
                }
            }

            if (needsMultiplication) tokens.push(new Token(TokenType.ASTERISK, "*"));
            tokens.push(nextToken);
        }

        this.mPunctuators = tokens.slice();
    }

    getTokens() {
        return this.mPunctuators;
    }

}