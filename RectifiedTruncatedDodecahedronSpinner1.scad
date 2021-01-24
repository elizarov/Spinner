$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 34; // [0:90]
// spinner height
sh = 30;
// base width
bw = 11;
// side width
sw = 3.5;
// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

g2 = gap / cos(a);

include <polyhedra.scad>;

spinner();

difference() {
    poly_wire(rectified(truncated_dodecahedron), sh, sw, fill_reg = [3]);
    spinner(g2, eps);
}


module spinner(dr = 0, eps = 0) {
    d = bw + 2 * dr;
    translate([0, 0, -eps]) {
        cylinder(d1 = d + 2 * eps, d2 = d + sh * tan(a), h = sh / 2 + eps);
        translate([0, 0, sh / 2])
            cylinder(d1 = d + sh * tan(a), d2 = d + 2 * eps, h = sh / 2 + eps);
    }
}

