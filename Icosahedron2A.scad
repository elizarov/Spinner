// shape height
sh = 30;
// plane width
pw = 1.5;
// wire width
bw = 3;

include <polyhedra.scad>;

translate([0, 0, sh / 2]) 
    face_rotate(icosahedron) 
        icosahedron2();

ivert = icosahedron[0];
iface = icosahedron[1];

ifr = face_dist(icosahedron);

module icosahedron0p() {
    poly_fill0(icosahedron, pw);
}

module icosahedron2() {
    ps = (sh - pw) / 2 / ifr;    
    module h4(a, b, c, d) {
        hull() {
            translate(ivert[a] * ps) icosahedron0p();
            translate(ivert[b] * ps) icosahedron0p();
            translate(ivert[c] * ps) icosahedron0p();
            translate(ivert[d] * ps) icosahedron0p();
        }
    }
    h4(0, 1, 10, 11);
    h4(4, 5, 6, 7);
    h4(2, 3, 8, 9);
    poly_wire0(icosahedron);
}
