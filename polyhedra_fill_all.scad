include <polyhedra.scad>;

// circumscribed diameter
sh = 20;
// number per row
m = 6;
// layout distance
ld = 2;

p = [
    // 5 Platonic Solids
    tetrahedron,
    octahedron,
    cube,
    icosahedron,
    dodecahedron,
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
    // 13 Catalan Solids
    triakis_tetrahedron,
    rhombic_dodecahedron,
    triakis_octahedron,
    tetrakis_hexahedron,
    deltoidal_icositetrahedron,
    disdyakis_dodecahedron,
    pentagonal_icositetrahedron,
    rhombic_triacontahedron,
    triakis_icosahedron,
    pentakis_dodecahedron,
    deltoidal_hexecontahedron,
    disdyakis_triacontahedron,
    pentagonal_hexecontahedron
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
        poly_fill(p[i], sh, circ=true);
}