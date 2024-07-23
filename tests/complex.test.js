

const { complex, Complex } = require("../math/complex.js");
const snapShotDiff = require("snapshot-diff");

expect.addSnapshotSerializer(snapShotDiff.getSnapshotDiffSerializer());


test("1 + 1 = 2", () => {
    expect(complex(1, 0).add(complex(1, 0)).re).toBe(2);
    expect(complex(1, 0).add(complex(1, 0)).im).toBe(0);
});


test("test complex addition", () => {
    expect(complex(1.12, 43.5).add(complex(51.309, 12.2))).toMatchSnapshot();
});