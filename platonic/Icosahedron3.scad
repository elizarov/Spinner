// shape height
sh = 30;
// border width
bw = 3;

include <../polyhedra.scad>;

poly_wire(icosahedron, fid = 1, fill_face_ids = [0, 4, 5, 14, 15, 19]);
