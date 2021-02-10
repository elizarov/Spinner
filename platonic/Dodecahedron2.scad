// shape heights
sh  = 26;
sh2 = 20;
// border width
bw = 3;

include <../polyhedra.scad>;

poly_wire(dodecahedron, sh);
translate([0, 0, (sh - sh2) / 2])
    poly_wire(snub_cube, sh2);    
