

const { complex, Complex } = require("../math/complex.js");


test("1 + 1 = 2", () => {
    expect(complex(1, 0).add(complex(1, 0)).re).toBe(2);
    expect(complex(1, 0).add(complex(1, 0)).im).toBe(0);
});