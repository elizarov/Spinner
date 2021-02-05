$fs = $preview ? 1 : 0.1;
$fa = $preview ? 10 : 2;

/* [Tolerance] */

// overhang angle
a = 40; // [0:90]
// gap between moving parts
gap = 0.45;
g2 = gap / cos(a); // gap for diagonal parts
// wall width 
ww = 1.6; // 4 x nozzle
// drilled diameter in the center (for strength)
dd = 0.6;
// flat surface height (drill hole below surface)
fsh = 0.6;

/* [Joint] */

// shape height
sh = 30;

// neck diameter
nd = dd + 2 * ww; 
echo("neck diameter", nd = nd);

// ball overlap length
bol = 0.5;
// ball diameter
bd = nd + 2 * gap + 2 * bol;
echo("ball diameter", bd = bd);

// outer neck diameter
ond = nd + 2 * gap; 
echo("outer neck diameter", ond = ond);

// neck housing diameter
nhd = ond + 2 * ww;
echo("neck housing diameter", nhd = nhd);

// ball housing diameter
bhd = bd + 2 * g2 + 2 * ww;
echo("ball housing diameter", bhd = bhd);

// surface chamfer length
scl = 1.5;

// outer surface diameter
osd = sh * tan(a) - (bhd + 2 * g2);
echo("outer surface diameter", osd = osd);

// real surface diameter
rsd = osd - 2 * scl;
echo("real surface diameter", rsd = rsd);
assert(rsd > nhd);

/* [Arm] */

// arm width
aw = 5;
assert(nhd > aw);
// arm length
al = 20;
// last arm length
lal = 15;

// overlap for clean differences and unions
eps = 1; 

/* [Layout] */
// number ofarms
na = 2;

for (i = [0:na - 1]) {
    last = i == na - 1;
    translate([i * al, 0, 0]) {
        arm(al = last ? lal : al, spindle = !last);
    }
}
#spindle();

module spindle() {
    difference() {
        spindle0();
        translate([0, 0, fsh])
            cylinder(d = dd, h = sh - 2 * fsh);
    }   
}

function mirror_y_2d(p) =
    let(m = [for(i = [len(p) - 2:-1:0]) [p[i].x, -p[i].y]])
        concat(p, m);

module spindle0(dr = 0, gap = 0, stand = true) {
    nx = nd / 2 + dr + gap;
    bx = bd / 2 + dr + gap / cos(a);
    sx = stand ? osd / 2 + dr + gap / cos(a) : nx;
    sy = sh / 2;
    bny = (bx - nx) / tan(a);
    sny = sy - (sx - nx) / tan(a);
    cx = stand ? sx - scl : nx;
    cy = sy - scl / tan(a);
    translate([0, 0, sh / 2])
        rotate_extrude()
            polygon(mirror_y_2d([
                [0, sy], [cx, sy], [cx, cy], [nx, sny], [nx, bny], 
                [bx, 0],
            ]));
}

module arm(al = al, housing = true, spindle = true) {
    xl = spindle ? al - nhd / 2 - gap : al;
    difference() {
        union() {
            translate([0, -aw / 2, 0])
                cube([xl, aw, sh]);
            if (housing) {
                spindle0(ww, gap, stand = false);
            }
        }
        if (housing) {
            wedge_bottom();    
            translate([0, 0, sh])
                rotate([180, 0, 0])
                    wedge_bottom();    
            spindle0(0, gap);
        }
        if (spindle) {
            translate([al, 0, 0]) 
                spindle0(ww, 2 * gap, stand = false);
        }
    }
    if (spindle) {
        translate([al, 0, 0]) 
            spindle();
    }
}

module wedge_bottom() {
    y = sh / 2;
    x = y * tan(a);
    rotate([90, 0, 0])
        translate([-bhd / 2, 0, -eps - bhd / 2])
            linear_extrude(bhd + 2 * eps)
                polygon([
                    [-eps, -eps], 
                    [x + eps * tan(a), -eps], 
                    [-eps, y + eps / tan(a)]
                ]);
}

module spinner(dr = 0, eps = 0, cdir = -1) {
    rotate_extrude() 
        polygon([]);
}
