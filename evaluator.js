/**
 * I think per-pixel function evaluation in javascript is fast enough to run in real time; rendering is the bottleneck
 * so JS generates the mesh and coloring can be done with a dynamic property passed as a uniform in reGL
 */





function evaluate(ast) {
    if (ast instanceof AssignExpression) {

    } else {
        
    }
}



function getEvaluatable(ast) {
    return function(z) {
        return z;
    }    
}