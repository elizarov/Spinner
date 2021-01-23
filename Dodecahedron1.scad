// shape height
sh = 30;
// border width
bw = 3;

translate([0, 0, sh / 2])
    dface_rotate()
        dodecahedron1();

// ------------ dodecahedron ------------

phi = (sqrt(5) + 1) / 2;

// dodecahedron vertices

dvert = [
    [-1/phi, 0, -phi], // 0
    [1/phi, 0, -phi], // 1
    [-1, -1, -1], // 2
    [-1, 1, -1], // 3
    [1, -1, -1], // 4
    [1, 1, -1], // 5
    [0, -phi, -1/phi], // 6
    [0, phi, -1/phi], // 7
    [-phi, -1/phi, 0], // 8
    [-phi, 1/phi, 0], // 9
    [phi, -1/phi, 0], // 10
    [phi, 1/phi, 0], // 11
    [0, -phi, 1/phi], // 12
    [0, phi, 1/phi], // 13
    [-1, -1, 1], // 14
    [-1, 1, 1], // 15
    [1, -1, 1], // 16
    [1, 1, 1], // 17
    [-1/phi, 0, phi], // 18
    [1/phi, 0, phi] // 19    
];

dface = [
    [0, 1, 5, 7, 3],
    [1, 0, 2, 6, 4],
    [0, 3, 9, 8, 2],
    [1, 4, 10, 11, 5],
    [8, 14, 12, 6, 2],
    [6, 12, 16, 10, 4],
    [7, 13, 15, 9, 3],
    [13, 7, 5, 11, 17],
    [8, 9, 15, 18, 14],
    [11, 10, 16, 19, 17],
    [18, 19, 16, 12, 14],
    [19, 18, 15, 13, 17]
];

df0c = (dvert[1] + dvert[0] + dvert[2] + dvert[6] + dvert[4]) / 5;
dfr = norm(df0c);
dfa = atan(df0c.y / df0c.z);

module dodecahedron1() {
    s = (sh - bw) / 2 / dfr;
    for (f = dface) {
        for (i = [0:4]) {
            j = (i + 1) % 5;
            p = f[i];
            q = f[j];
            if (p < q) {
                hull() {
                    translate(dvert[p] * s) dodecahedron0();
                    translate(dvert[q] * s) dodecahedron0();
                }                
            }
        }
    }   
}

module dodecahedron0() {
    scale(bw / 2 / dfr)
        polyhedron(dvert, dface);
}

module dface_rotate() {
    rotate([dfa, 0, 0]) children();
}