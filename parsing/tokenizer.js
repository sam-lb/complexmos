

const ALPHA = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
const NUM = "1234567890";
const OPERATORS = "+-*/^";
const ACCESSOR_OR_DEC = ".";
const OPEN_PAREN = "(";
const CLOSE_PAREN = ")";
const OPEN_BRACKET = "[";
const CLOSE_BRACKET = "]";
const ARG_AND_ITEM_SEP = ",";
const EQUALS = "=";
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
        argAndItemSep: 8,
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
    let numBuffer = "";
    let readingDecimalPart = false;
    const tokens = [];

    const clearIdentifierBuffer = (precedesImplicit=false) => {
        const bufEmpty = buffer.length === 0;
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
        if (!bufEmpty) {
            for (let identifier of identifiers) {
                tokens.push(new Token(identifier, Token.types.identifier));
                tokens.push(new Token("*", Token.types.operator));
            }
            if (!precedesImplicit) tokens.pop();
        }
        return true;
    };

    const clearNumberBuffer = (precedesImplicit=false) => {
        if (numBuffer.length > 0) {
            const value = parseFloat(numBuffer);
            if (isNaN(value)) {
                tracker.error(`Invalid number literal ${numBuffer}`);
                return false;
            }
            tokens.push(new Token(complex(value, 0), Token.types.number));
            if (precedesImplicit) {
                tokens.push(new Token("*", Token.types.operator));
            }
        }
        numBuffer = "";
        readingDecimalPart = false;
        return true;
    };

    const possiblyValidLeftOperand = () => {
        return (tokens.length > 0 && (tokens[tokens.length - 1] === CLOSE_BRACKET || tokens[tokens.length - 1] === CLOSE_PAREN));
    };

    // main tokenization loop
    for (let j=0; j<text.length; j++) {
        const character = text[j];
        if (ALPHA.includes(character)) {                                                // identifier character
            if (numBuffer.length > 0) {                                                 // beginning of identifier, end of number
                if (!clearNumberBuffer(true)) return null;                              // clear number buffer
            }
            buffer += character;                                                        // add character to identifier buffer
        } else if (NUM.includes(character)) {                                           // number character
            if (buffer.length > 0) {
                tracker.error("Implicit multiplication of the form {identifier}{number} (e.g. x2) is not supported.");
                return null;
            }
            numBuffer += character;
        } else if (character === ACCESSOR_OR_DEC) {
            if (numBuffer.length > 0) {                                                 // decimal point
                if (readingDecimalPart) {
                    tracker.error("Multiple decimal points in number");
                    return null;
                }
                numBuffer += character;
                readingDecimalPart = true;
            } else if (buffer.length > 0) {                                             // accessor

            } else {
                tracker.error(`Unexpected token ${ACCESSOR_OR_DEC}`);
                return null;
            }
        } else if (character === "-") {
            if (buffer.length > 0) {                                                    // subtraction operator
                if (!clearIdentifierBuffer()) return null;
                tokens.push(new Token("-", Token.types.operator));
            } else if (numBuffer.length > 0) {
                if (!clearNumberBuffer()) return null;
                tokens.push(new Token("-", Token.types.operator));
            } else if (possiblyValidLeftOperand()) {
                tokens.push(new Token("-", Token.types.operator));
            } else {                                                                    // unary negation operator
                tokens.push(new Token(complex(-1, 0), Token.types.number));
                tokens.push(new Token("*", Token.types.operator));
            }
        } else if (OPERATORS.includes(character)) {
            if (buffer.length > 0 || numBuffer.length > 0 || possiblyValidLeftOperand()) {
                if (!clearIdentifierBuffer()) return null;
                if (!clearNumberBuffer()) return null;
            } else {
                tracker.error(`Unexpected operator ${character}`);
                return null;
            }
            tokens.push(new Token(character, Token.types.operator));
        } else if (character === ARG_AND_ITEM_SEP) {

        } else if (character === OPEN_PAREN) {

        } else if (character === CLOSE_PAREN) {
            
        } else if (character === OPEN_BRACKET) {

        } else if (character === CLOSE_BRACKET) {

        } else {
            tracker.error(`Unexpected token ${character}`);
            return null;
        }
    }
    if (!clearIdentifierBuffer()) return null;
    if (!clearNumberBuffer()) return null;
    
    return {
        tokens: tokens,
        expressionType: null,
    };
}