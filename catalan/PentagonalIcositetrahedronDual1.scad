// shape height
sh = 50;
// border width
bw = 3;
// connector width
cw = 2;

include <../polyhedra.scad>;

poly_wire_dual(pentagonal_icositetrahedron, dual_fill_face_sizes = [4]);
