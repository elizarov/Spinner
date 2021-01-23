
// shape height
sh = 30;
// border width
bw = 3;


translate([0, 0, sh / 2]) 
    oface_rotate() 
        octahedron1();



// ------------ octahedron ------------

// octahedron vertices
overt = [
    [0, 0, -1], // 0
    [1, 0, 0], // 1
    [0, 1, 0], // 2
    [-1, 0, 0], // 3
    [0, -1, 0], // 4
    [0, 0, 1] // 5
];

oface = [
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 4],
    [0, 4, 1],
    [5, 2, 1],
    [5, 3, 2],
    [5, 4, 3],
    [5, 1, 4]
];

of0c = (overt[0] + overt[1] + overt[2]) / 3;
ofr = norm(of0c);
ofax = atan(of0c.y / of0c.z);

of0rx = of0c * [
    [1, 0, 0],
    [0, cos(ofax), -sin(ofax)],
    [0, sin(ofax), cos(ofax)]
];


ofay = atan(of0rx.x / of0rx.y);

module oface_rotate() {    
    rotate([ofax, ofay, 0]) children();
}

module octahedron0() {
    scale(bw / 2 / ofr) polyhedron(overt, oface);    
}

module octahedron1() {
    s = (sh - bw) / 2 / ofr;
    for (f = oface) {
        for (i = [0:2]) {
            j = (i + 1) % 3;
            p = f[i];
            q = f[j];
            if (p < q) {
                hull() {
                    translate(overt[p] * s) octahedron0();
                    translate(overt[q] * s) octahedron0();
                }                
            }
        }
    }        
}