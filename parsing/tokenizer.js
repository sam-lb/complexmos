

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



const expressionTypes = {
    constantAssignment: 0,
    functionAssignment: 1,
    literal: 2,
    expression: 3,
};

const dataTypes = {
    number: 0,
    matrix: 1,
    function: 2,
    array: 3,
};


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

    toString() {
        return `{${this.data}}`
    }

}


function buildSearchTrieFromScope(scope) {
    return new Trie(Object.keys(scope["builtin"]));
}


function tokenize(text, tracker, scope) {
    console.log(`tokenizing ${text}...`);
    if (text.length === 0) {
        tracker.error("Cannot tokenize empty string");
        return null;
    }

    if (text.includes(RESERVED)) {
        tracker.error(`text contains reserved character ${RESERVED}`);
        return null;
    }

    const identifierLookup = buildSearchTrieFromScope(scope);
    let buffer = "";
    let numBuffer = "";
    let readingDecimalPart = false;
    const tokens = [];
    let expectFunctionCall = false;
    let lastToken = null;

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
        for (let i=0; i<identifiers.length; i++) {
            const isFunction = scope.builtin[identifiers[i]].dataType === dataTypes.function;
            if (i === identifiers.length - 1) {
                expectFunctionCall = isFunction;
            } else {
                tracker.error("Function calls require parentheses.");
                return false;
            }
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
        return (
            tokens.length > 0 && (tokens[tokens.length - 1].type === Token.types.closeBracket
                                || tokens[tokens.length - 1].type === Token.types.closeParen)
        );
    };

    const requiresImplicit = () => {
        if (tokens.length === 0) return false;
        const token = tokens[tokens.length-1];
        return (
            token.type === Token.types.closeBracket 
            || token.type === Token.types.closeParen
            || token.type === Token.types.identifier
            || token.type === Token.types.number
        );
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
            if (!clearIdentifierBuffer()) return null;
            if (!clearNumberBuffer()) return null;
            if (requiresImplicit() && !expectFunctionCall) {
                tokens.push(new Token("*", Token.types.operator));
            }
            tokens.push(new Token(OPEN_PAREN, Token.types.openParen));
            expectFunctionCall = false;
        } else if (character === CLOSE_PAREN) {
            if (!clearIdentifierBuffer()) return null;
            if (!clearNumberBuffer()) return null;
            tokens.push(new Token(CLOSE_PAREN, Token.types.closeParen));
        } else if (character === OPEN_BRACKET) {
            if (!clearIdentifierBuffer()) return null;
            if (!clearNumberBuffer()) return null;
            if (requiresImplicit()) {
                tokens.push(new Token("*", Token.types.operator));
            }
            tokens.push(new Token(OPEN_BRACKET, Token.types.openBracket));
        } else if (character === CLOSE_BRACKET) {
            if (!clearIdentifierBuffer()) return null;
            if (!clearNumberBuffer()) return null;
            tokens.push(new Token(CLOSE_BRACKET, Token.types.CLOSE_BRACKET));
        } else {
            tracker.error(`Unexpected token ${character}`);
            return null;
        }
        if (expectFunctionCall) {
            tracker.error(`Function calls require parentheses`);
            return null;
        }
        lastToken = tokens[tokens.length-1];
    }
    if (!clearIdentifierBuffer()) return null;
    if (!clearNumberBuffer()) return null;
    if (expectFunctionCall) {
        tracker.error(`Function calls require parentheses`);
        return null;
    }
    
    return {
        tokens: tokens,
        expressionType: null,
    };
}