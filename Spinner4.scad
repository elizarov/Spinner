

$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 40; // [0:90]
// spinner height
sh = 30;
// base width
bw = 10;

// spinner shell width
sw = 6;
// spinner shell length
sl = 21;

// vertical connector width
vw = 3.5; 
// verical connector inner radium
vr = 12; 

ch = vw / (2 + sqrt(2)); // chamfer for vertical connectors

// support pillar radius
spr = 16;
// support pillar height
sph = 13;

// number of pellets
n = 2;
// number of pieces in vertical connectors
k = 3;

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

vc();


module shell() {
    for (i = [0:n-1])
        rotate([0, 0, i * 360 / n])
            translate([0, -sw / 2, 0])
                chamfer_cube([sl, sw, sh]);
    vc();
    translate([0, 0, sh])
        rotate([180, 0, 0]) vc();
}


module vc() {
    vl = 2 * (vr + vw / 2) * tan(180 / k / n);
    b = 360 / k / n;
    for (i = [0:k * n - 1])
        rotate([0, 0, (i + 0.5) * 360 / k / n]) {
            intersection() {
                translate([vr - vw / 2, - vl / 2, 0])
                    chamfer_cube_2([vw, vl, vw]);
                sector(b, vr + vw, sh);
            }
            if (i % k == 0 || i % k == k - 1) 
                connection(b, i % k == 0 ? 1 : -1);
        }
}

function vnorm(v) = v / norm(v);

function rotate_matrix(a0, b0) = 
    let(
        a = vnorm(a0),
        b = vnorm(b0),
        v = cross(a, b),
        s = norm(v),
        c = a * b,
        vm = [
            [0, -v.z, v.y],
            [v.z, 0, -v.x],
            [-v.y, v.x, 0]
        ],
        im = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ],
        rm = im + vm + vm * vm * ((1 - c) / (s * s))
    )
        [concat(rm.x, 0), concat(rm.y, 0), concat(rm.z, 0)];

module connection(b, d) {
    x1 = spr;
    y1 = x1 * tan(b / 2);
    z1 = sph;
    x2 = vr;
    y2 = x2 * tan(b / 2);
    z2 = vw / 2;
    p1 = [x1, -d * y1, z1];
    p2 = [x2, d * y2, z2];
    cl = norm(p1 - p2);
    translate(p1)
        multmatrix(rotate_matrix([0, 1, 0], p2 - p1))
            translate([-vw / 2, 0, -vw / 2])
                chamfer_cube_2([vw, cl, vw]);
}

module sector(b, r, h) {
    y = r * tan(b / 2);
    linear_extrude(h) {
        polygon([[0, 0], [r, y], [r, -y]]);
    }
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

module chamfer_cube_2(s, ch = ch) {
    difference() {
        cube(s);
        chamfer_cube_bottom_2(s, ch);
        translate([0, s.y, s.z])
            rotate([180, 0, 0]) chamfer_cube_bottom_2(s, ch);
    }
}

module chamfer_cube_bottom_2(s, ch) {
    translate([s.x, 0, 0])
        rotate([0, 0, 90]) chamfer(s.y, ch);
    translate([0, s.y, 0])
        rotate([0, 0, -90]) chamfer(s.y, ch);
}


module chamfer_cube_bottom(s, ch) {
    chamfer_cube_bottom_2(s, ch);
    chamfer(s.x, ch);
    translate([s.x, s.y, 0])
        rotate([0, 0, 180]) chamfer(s.x, ch);
    translate([0, 0, s.z])
        rotate([0, 90, 0]) chamfer(s.z, ch);
    translate([s.x, 0, 0])
        rotate([0, -90, 0]) chamfer(s.z, ch);
}

