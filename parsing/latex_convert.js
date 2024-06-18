function clearFractions(text, depth=0) {
    /*
    convert \frac{num}{denom} to ((num)/(denom))
    */
    if (depth > 9) {
        console.error("fractions nested too deeply");
        return "";
    }
    let newString = "", letter, buffer, bCount, j, fracOpen=false;
    for (let i=0; i<text.length; i++) {
        letter = text[i];
        if (i >= 4 && text.slice(i-4, i) === "frac") fracOpen = true;
        if (letter === "{" && fracOpen) {
            buffer = "";
            bCount = 1;
            for (j=i+1; j<text.length; j++) {
                if (text[j] === "{") {
                    bCount++;
                } else if (text[j] === "}") {
                    bCount--;
                if (bCount === 0) break;
                }
                buffer += text[j];
            }
            if (i >= 4 && text.slice(i-4, i) === "frac") {
                newString += "((" + clearFractions(buffer, depth+1) +")/"
            } else {
                newString += "(" + clearFractions(buffer, depth+1) + "))";
                fracOpen = false;
            }
            i = j;
        } else {
            newString += letter;
        }
    }
    return newString.replace(/\\frac/g, "");
}

/*
sum and prod will be written as functions e.g.
sum(k, [1,...,infinity], 1/k^2)
prod(p, [2, 3, 5, 7, 11, 13], 1/(1-1/p))
*/
function cleanLatex(text) {
    text = text.replace(/\s/g, "");                                                                 // remove whitespace
    text = text.replace(/\\cdot/g, "*").replace(/\\operatorname{(.*?)}/g, (tot, group1)=>group1);   // remove operator wrappers
    text = clearFractions(text);                                                                    // convert latex fractions to a/b
    text = text.replace(/(\\left)|(\\right)/g, "");                                                 // remove \left and \right
    text = text.replace(/\\/g, "");                                                                 // remove backslashes
    text = text.replace(/{/g, "(").replace(/}/g, ")");                                              // replace curly brackets with parentheses
    return text;
}


/*
test strings

4\cos\left(a+bi\right)\cdot\sin\left(x\right)+\frac{1}{2iz}
\frac{1}{1+\frac{1}{x}+\frac{1}{y}}
*/



module.exports = {
    cleanLatex,
};