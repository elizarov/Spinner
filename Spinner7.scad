$fs = $preview ? 1 : 0.2;
$fa = $preview ? 10 : 2;

include <bezier.scad>;

// angle
a = 45; // [0:90]
// spinner height
sh = 10;
// neck width
nw = 6;
// base width
bw = nw + sh * sin(a);
echo("base width", bw = bw);

// spinner shell diameter
ds = bw + 9;
rs = ds / 2; // radius

// arm radius
ra = 35; 
// arm width
wa = 12;

// number of arms
n = 3;

// shell chamfers
ch = 1.5;
// spinner chamfers
sch = 1.5;

echo("real base width", bw - 2 * sch * tan(a));

// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

// drilled diameter in the center (for strength)
dd = 0.6;
// flat surface height (drill hole below surface)
fsh = 0.6;

g2 = gap / cos(a);

// compute arm shape
b = 180 / n;
x0 = rs * cos(b);
y0 = rs * sin(b);
y1 = wa / 2;
x1 = x0 + (y0 - y1) * tan(b);
cn = [cos(b), sin(b)];

function arm_slice(x, y, n) = 
    let(cn = ch * n)
        [
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
    spinner0(g2, eps);
    for (i = [0:n-1])
        rotate([0, 0, 360 * i / n])
            translate([ra, 0, 0]) 
               spinner0(g2, eps);
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

module spinner() {
    difference() {
        spinner0(cy = sch);
        translate([0, 0, fsh])
            cylinder(d = dd, h = sh - 2 * fsh);
    }
}

module spinner0(dr = 0, eps = 0, cy = 0) {
    ro = bw / 2 + dr;
    rn = nw / 2 + dr;
    cx = cy * tan(a);
    rotate_extrude() 
        polygon([
            [ro - cx, -eps],
            [ro - cx, cy],
            [rn, sh / 2],
            [ro - cx, sh - cy],
            [ro - cx, sh + eps],
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

