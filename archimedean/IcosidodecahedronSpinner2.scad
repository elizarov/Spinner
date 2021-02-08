$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 35; // [0:90]
// spinner height
sh = 30;
// base width
bw = 10;
// side width
sw = 3.5;

// vertical chamfer displacement
cv = 1.3;
// bottom chamfer displacement (total)
cb = 1.9;

// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1;

g2 = gap / cos(a);

include <../polyhedra.scad>;

spinner();

difference() {
    poly_wire(icosidodecahedron, sh, sw, fill_face_sizes = [3]);
    spinner(g2, eps);
    chamfer_both(bw + 2 * g2);
}


module spinner(dr = 0, eps = 0) {
    d = bw + 2 * dr;
    translate([0, 0, -eps]) {
        cylinder(d1 = d + 2 * eps, d2 = d + sh * tan(a), h = sh / 2 + eps);
        translate([0, 0, sh / 2])
            cylinder(d1 = d + sh * tan(a), d2 = d + 2 * eps, h = sh / 2 + eps);
    }
}

module chamfer_both(d) {
    chamfer_bottom(d);
    translate([0, 0, sh])
        rotate([180, 0, 0])
            chamfer_bottom(d);
}

module chamfer_bottom(d) {
    translate([0, 0, -eps])
        cylinder(d1 = d + 2 * cb + 2 * eps, d2 = d, h = cb + eps);
    cylinder(d = d + 2 * cv, h = cb + eps);
}
