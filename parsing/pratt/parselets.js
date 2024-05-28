


class InfixParselet {

    parse(parser, left, token) {

    }

    getPrecedence() {

    }

}


class PrefixParselet {

    parse(parser, token) {

    }

    getPrecedence() {

    }

}


class AssignParselet extends InfixParselet {

    parse(parser, left, token) {
        const right = parser.parseExpression(Precedence.ASSIGNMENT - 1);
        if (!(left instanceof NameExpression || left instanceof CallExpression)) {
            console.error("lhs of assignment should be identifier or function name with arguments");
        }

        return new AssignExpression(left.toString(), right);
    }

    getPrecedence() {
        return Precedence.ASSIGNMENT;
    }

}


class BinaryOperatorParselet extends InfixParselet {

    constructor(precedence, mIsRight) {
        super();
        this.mPrecedence = precedence;
        this.mIsRight = mIsRight;
    }

    parse(parser, left, token) {
        const right = parser.parseExpression(this.mPrecedence - (this.mIsRight ? 1 : 0));
        return new OperatorExpression(left, token.mtype, right);
    }

    getPrecedence() {
        return this.mPrecedence;
    }

}


class CallParselet extends InfixParselet {

    constructor() {
        super();
    }

    parse(parser, left, token) {
        const args = [];
        
        if (!parser.peek(TokenType.RIGHT_PAREN)) {
            while (true) {
                args.push(parser.parseExpression());
                if (parser.peek(TokenType.RIGHT_PAREN)) {
                    parser.consume(TokenType.RIGHT_PAREN);
                    break;
                }
                parser.consume(TokenType.COMMA);
            }
        }

        return new CallExpression(left, args);
    }

    getPrecedence() {
        return Precedence.CALL;
    }

}


class GroupParselet extends PrefixParselet {

    parse(parser, token) {
        const expr = parser.parseExpression();
        parser.consume(TokenType.RIGHT_PAREN);
        return expr;
    }

}


class NameParselet extends PrefixParselet {

    constructor() {
        super();
    }

    parse(parser, token) {
        return new NameExpression(token.text);
    }

}


class PostfixOperatorParselet extends InfixParselet {

    constructor(precedence) {
        super();
        this.mPrecedence = precedence;
    }

    parse(parser, left, token) {
        return new PostfixExpression(token.mtype, left);
    }

    getPrecedence() {
        return this.mPrecedence;
    }

}


class PrefixOperatorParselet extends PrefixParselet {

    constructor(precedence) {
        super();
        this.mPrecedence = precedence;
    }

    parse(parser, token) {
        const right = parser.parseExpression(this.mPrecedence);
        return new PrefixExpression(token.mtype, right);
    }

}