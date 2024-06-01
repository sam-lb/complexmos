

class ExpressionParser extends Parser {

    constructor(tokens) {
        super(tokens);

        this.registerPrefix(TokenType.NUMBER, new NumberParselet());
        this.registerPrefix(TokenType.NAME, new NameParselet());

        this.registerInfix(TokenType.ASSIGN, new AssignParselet());
        this.registerPrefix(TokenType.LEFT_PAREN, new GroupParselet());
        this.registerInfix(TokenType.LEFT_PAREN, new CallParselet());

        this.prefix(TokenType.PLUS, Precedence.PREFIX);
        this.prefix(TokenType.MINUS, Precedence.PREFIX);

        this.infixLeft(TokenType.COLON, Precedence.CALL);
        this.infixLeft(TokenType.PLUS, Precedence.SUM);
        this.infixLeft(TokenType.MINUS, Precedence.SUM);
        this.infixLeft(TokenType.ASTERISK, Precedence.PRODUCT);
        this.infixLeft(TokenType.SLASH, Precedence.PRODUCT);
        this.infixRight(TokenType.CARET, Precedence.EXPONENT);

    }

    prefix(tokenType, precedence) {
        this.registerPrefix(tokenType, new PrefixOperatorParselet(precedence));
    }

    infixLeft(tokenType, precedence) {
        this.registerInfix(tokenType, new BinaryOperatorParselet(precedence, false));
    }

    infixRight(tokenType, precedence) {
        this.registerInfix(tokenType, new BinaryOperatorParselet(precedence, true));
    }

}