// shape height
sh = 50;
// border width
bw = 2.5;
// connector width
cw = 1.6;

include <../polyhedra.scad>;

poly_wire_dual(pentagonal_icositetrahedron, dual_fill_reg = [4]);
