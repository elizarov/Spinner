include <../polyhedra.scad>;

// circumscribed diameter
sh = 30;
// border width
bw = 3;
// number per row
m = 4;
// layout distance
ld = 2;

p = [
    // 13 Arhimedean Solids
    truncated_tetrahedron,
    cuboctahedron,
    truncated_cube,
    truncated_octahedron,
    rhombicuboctahedron,
    rhombitruncated_cuboctahedron,
    snub_cube,
    icosidodecahedron,
    truncated_dodecahedron,
    truncated_icosahedron,
    rhombicosidodecahedron,
    rhombitruncated_icosidodecahedron,
    snub_dodecahedron,
];

n = len(p);
last_row = floor((n - 1) / m);
first_row_c = floor((m + n - last_row * m) / 2);
skip1 = m - first_row_c;
last_row_c = n - last_row * m + skip1;
for(i = [0:n - 1]) {
    j = i + skip1;
    row = floor(j / m);
    ofs = 
        row == 0 ? -skip1 / 2 :
        row == last_row ? (m - last_row_c) / 2 :
        0;
    y = row * (sh + ld);
    x = (ofs + j % m) * (sh + ld);
    translate([x, y, 0])
        poly_wire(p[i], sh, bw, circ=true);
}