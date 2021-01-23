
// shape height
sh = 30;
// border width
bw = 3;


translate([0, 0, sh / 2]) 
    iface_rotate() 
        icosahedron1();



// ------------ icosahedron ------------

phi = (sqrt(5) + 1) / 2;

// icosahedron vertices
ivert = [
    [0, -1, -phi], // 0
    [0, 1, -phi], // 1
    [-phi, 0, -1], // 2
    [phi, 0, -1], // 3
    [-1, -phi, 0], // 4
    [-1, phi, 0], // 5
    [1, -phi, 0], // 6
    [1, phi, 0], // 7
    [-phi, 0, 1], // 8
    [phi, 0, 1], // 9
    [0, -1, phi], // 10
    [0, 1, phi] // 11    
];

iface = [
    [0, 1, 2],
    [1, 0, 3],
    [0, 2, 4],
    [2, 1, 5],
    [1, 3, 7],
    [3, 0, 6],
    [1, 7, 5],
    [0, 4, 6],
    [2, 8, 4],
    [2, 5, 8],
    [3, 6, 9],
    [3, 9, 7],
    [4, 10, 6],
    [5, 7, 11],
    [8, 10, 4],
    [5, 11, 8],
    [9, 11, 7],
    [6, 10, 9],
    [8, 11, 10],
    [9, 10, 11]
];

if0c = (ivert[0] + ivert[1] + ivert[2]) / 3;
ifr = norm(if0c);
ifa = atan(if0c.x / if0c.z);


    

module iface_rotate() { 
    rotate([0, ifa, 0]) children();
}

module icosahedron0() {
    scale(bw / 2 / ifr) polyhedron(ivert, iface);    
}

module icosahedron1() {
    s = (sh - bw) / 2 / ifr;
    for (f = iface) {
        for (i = [0:2]) {
            j = (i + 1) % 3;
            p = f[i];
            q = f[j];
            if (p < q) {
                hull() {
                    translate(ivert[p] * s) icosahedron0();
                    translate(ivert[q] * s) icosahedron0();
                }                
            }
        }
    }    
    module hface(i) {
        f = iface[i];
        hull() {
            translate(ivert[f[0]] * s) icosahedron0();
            translate(ivert[f[1]] * s) icosahedron0();
            translate(ivert[f[2]] * s) icosahedron0();
        }
    }
    hface(0);
    hface(4);
    hface(5);
    hface(14);
    hface(15);
    hface(19);
}