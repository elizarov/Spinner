// shape height
sh = 50;
// border width
bw = 2.5;
// dual width
dw = 1.5;

cw = 0;

include <../polyhedra.scad>;

fid = filter_face_ids(truncated_icosahedron, [5])[0];

poly_wire_dual(
    truncated_icosahedron,
    dw = dw,
    fill_face_sizes = [5],
    fid = fid
);
