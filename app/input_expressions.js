


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

}


class FunctionDefinition extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        this.requirements = [];        
    }

    readTokens() {
        this.name = null;
        this.arguments = {
            arg1: "type1",
        };
        this.definition = null;
        // this.ast = parse(this.tokens);
    }

}


class VariableDefinition extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        this.requirements = [];
    }

    readTokens() {
        this.name = null;
        this.definition = null;
        this.sliderable = false;
        this.type = "type";
    }

}


class EvaluatableLine extends InputExpression {

    constructor(tokens, id) {
        super(tokens, id);
    }

    findRequirements() {
        this.requirements = [];
    }

    readTokens() {
        
    }

}


module.exports = {
    FunctionDefinition,
    VariableDefinition,
    EvaluatableLine,
};