

$fs = $preview ? 0.5 : 0.2;
$fa = $preview ? 10 : 2;

// angle
a = 40; // [0:90]
// spinner height
sh = 30;
// neck width
nw = 6;
// base width
bw = nw + sh * sin(a);

// spinner shell width
sw = 4;
// spinner shell length
sl = 25;

// neck cylinder width
ncw = 12;
// neck cylinder height
nch = 10;

// number of pellets
n = 6;

// cube chamfers
ch = 1;
// bottom chamfer displacement
cb = 3;

// gap between moving parts
gap = 0.5;
// overlap for clean differences and unions
eps = 0.1; 

g2 = gap / cos(a);

spinner();

#intersection() {
    difference() {
        shell();
        spinner(g2, eps, 1);
    }
}


module shell() {
    translate([0, 0, (sh - nch) / 2])
        cylinder(d = ncw, h = nch);
    for (i = [0:n-1])
        rotate([0, 0, i * 360 / n])
            translate([0, -sw / 2, 0])
                chamfer_cube([sl, sw, sh]);
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

