// shape height
sh = 30;
// border width
bw = 5;

include <polyhedra.scad>;

translate([0, 0, sh / 2]) 
    face_rotate(icosahedron) 
        icosahedron2();

ivert = icosahedron[0];
iface = icosahedron[1];

ifr = face_dist(icosahedron);

ipd = diameter(icosahedron);

module icosahedron0() {
    poly_fill0(icosahedron, ipd, bw);
}

module icosahedron2() {
    s = (sh - bw) / 2 / ifr;
    module h4(a, b, c, d) {
        hull() {
            translate(ivert[a] * s) icosahedron0();
            translate(ivert[b] * s) icosahedron0();
            translate(ivert[c] * s) icosahedron0();
            translate(ivert[d] * s) icosahedron0();
        }
    }
    h4(0, 1, 10, 11);
    h4(4, 5, 6, 7);
    h4(2, 3, 8, 9);
}
