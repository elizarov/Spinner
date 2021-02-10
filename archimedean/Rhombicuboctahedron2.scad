// shape heights
sh  = 26;
sh2 = 19;
// border width
bw = 3;

include <../polyhedra.scad>;

poly_wire(rhombicuboctahedron, sh);
translate([0, 0, (sh - sh2) / 2])
    poly_wire(snub_dodecahedron, sh2);

