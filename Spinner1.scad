

$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 40; // [0:90]
// spinner height
sh = 30;
// base width
bw = 10;
// side width
sw = 3;
// side depth
sd = 4;
// vertical chamfer displacement
cv = 2;
// total chamfer displacement
ch = 3;
// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

g2 = gap / cos(a);
odr = (sd + gap) / cos(a);
or = bw / 2 + sh / 2 * tan(a) + odr; //gap / cos(a) + sd;

swz = sw * cos(a);

spinner();

intersection() {
    difference() {
        spinner(odr);
        spinner(g2, eps);
        for (i = [0:3])
            rotate([0, 0, i * 90])
                side_cut();
        chamfer_both(bw + 2 * g2);
    }
    //cylinder(r = or, h = sh);
}

*side_cut();

module spinner(dr = 0, eps = 0) {
    d = bw + 2 * dr;
    translate([0, 0, -eps]) {
        cylinder(d1 = d + 2 * eps, d2 = d + sh * tan(a), h = sh / 2 + eps);
        translate([0, 0, sh / 2])
            cylinder(d1 = d + sh * tan(a), d2 = d + 2 * eps, h = sh / 2 + eps);
    }
}

module side_cut() {
//    translate([sw / 2, sw / 2, swz])
//        cube([or, or, sh - 2 * swz]);
    d = sw / 2;
    e = or;
    polyhedron(
        points = [
            [d, d, swz], // 0
            [e, e, swz], // 1
            [d, e, sh / 2], // 2
            [e, d, sh / 2], // 3
            [d, d, sh - swz], // 4
            [e, e, sh - swz]  // 5
        ], 
        faces = [
            [0, 3, 1],
            [0, 4, 3],
            [3, 5, 1],
            [4, 5, 3],
            
            [0, 1, 2],
            [0, 2, 4],
            [1, 5, 2],
            [4, 2, 5]         
        ]
    );
}


module chamfer_both(d) {
    chamfer_bottom(d);
    translate([0, 0, sh])
        rotate([180, 0, 0]) 
            chamfer_bottom(d);
}

module chamfer_bottom(d) {
    translate([0, 0, -eps])
        cylinder(d1 = d + 2 * ch + 2 * eps, d2 = d, h = ch + eps);
    cylinder(d = d + 2 * cv, h = ch + eps);
}

