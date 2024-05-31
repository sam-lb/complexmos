



const pColorMap = {
    0: "#f00",
    1: "#0f0",
    2: "#00f",
    3: "#ff0",
    4: "#0ff",
    5: "#f0f",
    6: "#f80",
    7: "#0f8",
    8: "#f08",
    9: "#8f0",
    10: "#80f",
};


function parenthesis_match(string) {
    let pCount = 0;
    let result = "";
    for (const char of string) {
        if (char === "(") {
            const color = pColorMap[Math.min(pCount, 10)];
            pCount++;

            result += `<span style="color: ${color};font-weight:1000;">(</span>`;
        } else if (char === ")") {
            pCount--;
            const color = pColorMap[Math.min(pCount, 10)];

            result += `<span style="color: ${color};font-weight:1000;">)</span>`;
        } else {
            result += char;
        }
    }
    return result;
}

