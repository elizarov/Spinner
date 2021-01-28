// shape height
sh = 40;
// border width
bw = 2;
// dual width
dw = 1.5;

cw = 0;

include <polyhedra.scad>;

fid = filter_face_ids(truncated_icosahedron, [5])[0];

poly_wire_dual(
    truncated_icosahedron, 
    dw = dw,
    fill_reg = [5],
    fid = fid
);