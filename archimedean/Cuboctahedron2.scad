// shape heights
sh  = 26;
sh2 = 20;
// border width
bw = 3;

include <../polyhedra.scad>;

poly_wire(cuboctahedron, sh);
translate([0, 0, (sh - sh2) / 2])
    poly_wire(rhombicosidodecahedron, sh2);