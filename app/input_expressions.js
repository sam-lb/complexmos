

const { tracker } = require("../parsing/errors.js");
const { TokenType } = require("../parsing/pratt/tokentype.js");
const { scope } = require("./scope.js");
const { ExpressionParser } = require("../parsing/pratt/expression_parser.js");
const { NumberExpression, PrefixExpression } = require("../parsing/pratt/expressions.js");



class PlottableType {
    static NUMBER = 0; // point
    static DOMAIN = 1; // something that can be shown as a domain coloring
    static PARAMETRIC = 2;
    static NONE = 3; // no way to plot, e.g. functions of more than one variable
}


class InputExpression {

    constructor(tokens, id) {
        this.tokens = tokens;
        this.id = id;
        this.findRequirements();
        this.readTokens();
    }

    findRequirements() {
        this.requirements = [];
    }

    readTokens() {
        
    }

    buildAST() {
        this.ast = (new ExpressionParser(this.tokens).parseExpression());
    }

}


class FunctionDefinition extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        this.requirements = [];
        this.locals = [];
        if (this.tokens[0]?.mtype !== TokenType.NAME) {
            tracker.error("Function definition does not begin with valid identifier");
            return;
        } else {
            this.name = this.tokens[0].text;
        }

        let assignmentEncountered = false;
        for (const token of this.tokens.slice(1, this.tokens.length)) {
            if (token.mtype === TokenType.ASSIGN) {
                assignmentEncountered = true;
                continue;   
            }

            if (token.mtype !== TokenType.NAME) continue;

            if (assignmentEncountered) {
                if (!scope.builtin[token.text] && !this.locals.includes(token.text)) {
                    this.requirements.push(token.text);
                }
            } else {
                if (!scope.builtin[token.text]?.isType) {
                    this.locals.push(token.text);
                }
            }
        }
    }

    readTokens() {
        this.arguments = {
            arg1: "type1",
        };
        this.definition = null;
    }

    buildAST() {
        super.buildAST();
        if (!tracker.hasError) this.ast = this.ast.mRight;
    }

}


class VariableDefinition extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        if (this.tokens[0]?.mtype !== TokenType.NAME) {
            tracker.error("Variable definition does not begin with valid identifier");
            return;
        } else {
            this.name = this.tokens[0].text;
        }

        this.requirements = this.tokens.slice(2, this.tokens.length).filter(token => {
            return token.mtype === TokenType.NAME && (!scope.builtin[token.text] || scope.builtin[token.text].isParameter);
        }).map(token => token.text);
    }

    readTokens() {
        this.definition = null;
        this.sliderable = false;
        this.type = "type";
    }

    buildAST() {
        super.buildAST();
        if (!tracker.hasError) this.ast = this.ast.mRight;
    }

    sliderBounds(fields) {
        // return the bounds of a slider for this variable, if applicable (if not, return null)
        const bounds = {};
        if (this.ast instanceof NumberExpression || (this.ast instanceof PrefixExpression && this.ast.mRight instanceof NumberExpression)) {
            // needs a slider
            const sliderName = fields[this.id]["sliderName"];
            if (sliderName === this.name) {
                bounds.min = fields[this.id]["savedBounds"].min;
                bounds.max = fields[this.id]["savedBounds"].max;
            } else {
                bounds.min = Math.floor(this.ast.mNumber);
                bounds.max = Math.floor(this.ast.mNumber + 1);
            }

            fields[this.id]["savedBounds"] = bounds;
            fields[this.id]["sliderName"] = this.name;
            return bounds;
        } else {
            return null;
        }
    }

}


class EvaluatableLine extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        this.requirements = this.tokens.filter(token => {
            return token.mtype === TokenType.NAME && (!scope.builtin[token.text] || scope.builtin[token.text].isParameter);
        }).map(token => token.text);
    }

    readTokens() {
        this.type = null;
        this.plottableType = null;
    }

}

module.exports = {
    FunctionDefinition,
    VariableDefinition,
    EvaluatableLine,
};