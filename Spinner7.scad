$fs = $preview ? 1 : 0.2;
$fa = $preview ? 10 : 2;

include <bezier.scad>;

// angle
a = 30; // [0:90]
// spinner height
sh = 10;
// neck width
nw = 4;
// base shell width
bw = nw + sh * sin(a);
// spinner shell diameter
ds = bw + 7;
rs = ds / 2; // radius

// arm radius
ra = 30; 
// arm width
wa = 7;

// number of arms
n = 3;

// cube chamfers
ch = 1;

// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

g2 = gap / cos(a);

// compute arm shape
b = 180 / n;
x0 = rs * cos(b);
y0 = rs * sin(b);
y1 = wa / 2;
x1 = x0 + (y0 - y1) * tan(b);
cn = [cos(b), sin(b)];

function arm_slice(x, y, cn) = [
    [x - cn.x, y - cn.y, 0],
    [x, y, ch],
    [x, y, sh - ch],
    [x - cn.x, y - cn.y, sh],
    [x - cn.x, -y + cn.y, sh],
    [x, -y, sh - ch],
    [x, -y, ch],
    [x - cn.x, -y + cn.y, 0]
];

module arm() {
    bezier_surface([
        arm_slice(x0, y0, cn),
        arm_slice(x1, y1, [0, 1]),
        arm_slice(ra / 2, y1, [0, 1]),
        arm_slice(ra - x1, y1, [0, 1]),
        arm_slice(ra - x0, y0, [-cn.x, cn.y])
    ]);
}

difference() {
    full_shell();
    spinner(g2, eps, 1);
    for (i = [0:n-1])
        rotate([0, 0, 360 * i / n])
            translate([ra, 0, 0]) 
               spinner(g2, eps, 1);
}

spinner();
for (i = [0:n-1])
    rotate([0, 0, 360 * i / n])
        translate([ra, 0, 0]) 
            spinner();

module full_shell() {
    shell();
    for (i = [0:n-1])
        rotate([0, 0, 360 * i / n]) {
            translate([ra, 0, 0]) 
                shell();
            arm();
        }
}

module shell() {
    rotate_extrude() 
        polygon([
            [0, 0],
            [rs - ch, 0],
            [rs, ch],
            [rs, sh - ch],
            [rs - ch, sh],
            [0, sh]
        ]);
}


module spinner(dr = 0, eps = 0, cdir = -1) {
    ro = bw / 2 + dr;
    rn = nw / 2 + dr;
    rotate_extrude() 
        polygon([
            [ro + cdir * ch, -eps],
            [ro, ch],
            [rn, sh / 2],
            [ro, sh - ch],
            [ro + cdir * ch, sh + eps],
            [0, sh + eps],
            [0, -eps],            
        ]);
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

