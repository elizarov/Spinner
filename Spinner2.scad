

$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 40; // [0:90]
// spinner height
sh = 30;
// base width
bw = 10;

// spinner shell width
sw = 15;
// spinner shell length
sl = 40;

// cube chamfers
ch = 1;
// vertical chamfer displacement
cv = 2;
// bottom chamfer displacement (total)
cb = 3;

// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

g2 = gap / cos(a);

swz = sw * cos(a);

spinner();

intersection() {
    difference() {
        shell();
        spinner(g2, eps);
        chamfer_both(bw + 2 * g2);
    }
}

module shell() {
    shell0();
    rotate([0, 0, 90]) shell0();   
}

module shell0() {
    translate([-sl / 2, -sw / 2, 0])
        chamfer_cube([sl, sw, sh]);
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

module chamfer(d, ch = ch) {
    translate([-eps, 0, 0])
        rotate([90, 0, 0]) rotate([0, 90, 0])
            linear_extrude(d + 2 * eps)
                polygon([[-eps, -eps], [ch + eps, -eps], [-eps, ch + eps]]);
}

module chamfer_cube(s, ch = ch) {
    difference() {
        cube(s);
        chamfer_cube_bottom(s, ch);
        translate([0, s.y, s.z])
            rotate([180, 0, 0]) chamfer_cube_bottom(s, ch);
    }
}

module chamfer_cube_bottom(s, ch) {
    chamfer(s.x, ch);
    translate([s.x, s.y, 0])
        rotate([0, 0, 180]) chamfer(s.x, ch);
    translate([s.x, 0, 0])
        rotate([0, 0, 90]) chamfer(s.y, ch);
    translate([0, s.y, 0])
        rotate([0, 0, -90]) chamfer(s.y, ch);
    translate([0, 0, s.z])
        rotate([0, 90, 0]) chamfer(s.z, ch);
    translate([s.x, 0, 0])
        rotate([0, -90, 0]) chamfer(s.z, ch);
}

