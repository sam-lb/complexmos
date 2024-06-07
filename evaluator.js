/**
 * I think per-pixel function evaluation in javascript is fast enough to run in real time; rendering is the bottleneck
 * so JS generates the mesh and coloring can be done with a dynamic property passed as a uniform in reGL
 */





function evaluate(ast) {
    // at this point, the ast has gone through the lexer and parser's checks already
    if (ast instanceof AssignExpression) {
        // the expression is trying to define a function
        const left = ast.mLeft;
        let varName;
        if (left instanceof CallExpression) {
            // the expression is defining a function
            varName = left.mFunction.mName;
        } else {
            // the expression is defining a variable
            varName = left.mName;            
        }

        valueScope[varName] = getEvaluatable(ast.mRight, left.mArgs);
        return null;
    } else {
        return getEvaluatable(ast);
    }
}



function getEvaluatable(ast, args=null) {
    console.log(args);
    return function(z) {
        return z;
    }    
}