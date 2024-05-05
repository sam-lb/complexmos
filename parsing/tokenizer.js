

const ALPHA = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
const NUM = "1234567890";
const OPERATORS = "+-*/^";
const RESERVED = "@";


class Token {

    static types = {
        identifier: 0,
        number: 1,
        operator: 2,
        openParen: 3,
        closeParen: 4,
        openBracket: 5,
        closeBracket: 6,
        equals: 7,
    };

    constructor(data, type) {
        this.data = data;
        this.type = type;
    }

}


function tokenize(text, tracker, scope) {
    if (text.length === 0) {
        tracker.error("Cannot tokenize empty string");
        return null;
    }

    if (text.includes(RESERVED)) {
        tracker.error(`text contains reserved character ${RESERVED}`);
        return null;
    }

    const identifierLookup = new Trie(scope["builtin"]);
    let buffer = "";
    const tokens = [];

    const clearIdentifierBuffer = () => {
        const identifiers = [];
        while (buffer.length > 0) {
            let matchFound = false;
            let possibleIdentifier, i;
            for (i=buffer.length-1; i>=0; i--) {
                possibleIdentifier = buffer.slice(0, i+1);
                if (identifierLookup.containsKey(possibleIdentifier)) {
                    matchFound = true;
                    break;
                }
            }
            if (matchFound) {
                identifiers.push(possibleIdentifier);
                buffer = buffer.slice(i+1, buffer.length);
            } else {
                tracker.error(`Undefined identifier ${buffer}`);
                return false;
            }
        }
        for (let identifier of identifiers) {
            tokens.push(new Token(identifier, Token.types.identifier));
            tokens.push(new Token("*", Token.types.operator));
        }
        tokens.pop();
        return true;
    }

    const clearNumberBuffer = () => {

    }

    // main tokenization loop
    for (let j=0; j<text.length; j++) {
        const character = text[j];
        if (ALPHA.includes(character)) {
            buffer += character;
        } else {
            if (!clearIdentifierBuffer()) return null;
        }
    }
    if (!clearIdentifierBuffer()) return null;
    
    return tokens;
}