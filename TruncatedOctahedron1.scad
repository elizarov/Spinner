
// shape height
sh = 30;
// border width
bw = 3;

tr = 1/3; // truncation ratio

translate([0, 0, sh / 2]) 
    toface_rotate() 
        trucated_octahedron1();

// ------------ octahedron ------------

// octahedron vertices
overt = [
    [0, 0, -1], // 0
    [1, 0, 0], // 1
    [0, 1, 0], // 2
    [-1, 0, 0], // 3
    [0, -1, 0], // 4
    [0, 0, 1] // 5
];

oface = [
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 4],
    [0, 4, 1],
    [5, 2, 1],
    [5, 3, 2],
    [5, 4, 3],
    [5, 1, 4]
];

// ------------ truncated octahedron ------------

// vertice pairs
tovp = [
    for (f = oface) 
        for (i = [0:2])
            let(j = (i + 1) % 3)
                [f[i], f[j]]
];

tovert = [
    for (p = tovp) 
        let(a = overt[p[0]])
        let(b = overt[p[1]])
            a + tr * (b - a)
];
        
function tovn(u, v) =
    search([[u, v]], tovp)[0];         
        
// hexagonal faces        
toface_h = [
    for (f = oface) [
        for (t = [0:5])
            let(i = t / 2)
            let(j = (i + 1) % 3)
                t % 2 == 0 ? tovn(f[i], f[j]) : tovn(f[j], f[i])
    ]
];        
  
// rectangular faces (TBD -- not computed now)       
toface_r = [
/*        
    for (i = [0:7]) [
        for (f = oface) 
            for (j = [0:2])
                if (f[j] == i)
                    let(k = (j + 1) % 3)
                        tovn(i, f[k])
    ]
*/        
];        
        
toface = concat(toface_h, toface_r);
        

tof0c = (tovert[0] + tovert[1] + tovert[2]) / 3;
tofr = norm(tof0c);
tofax = atan(tof0c.y / tof0c.z);

tof0rx = tof0c * [
    [1, 0, 0],
    [0, cos(tofax), -sin(tofax)],
    [0, sin(tofax), cos(tofax)]
];


tofay = atan(tof0rx.x / tof0rx.y);

module toface_rotate() {    
    rotate([tofax, tofay, 0]) children();
}

module trucated_octahedron0() {
    scale(bw / 2 / tofr) hull() polyhedron(tovert, toface);    
}

module trucated_octahedron1() {
    s = (sh - bw) / 2 / tofr;
    for (f = toface_h) {
        for (i = [0:5]) {
            j = (i + 1) % 6;
            p = f[i];
            q = f[j];
            hull() {
                translate(tovert[p] * s) trucated_octahedron0();
                translate(tovert[q] * s) trucated_octahedron0();
            }                
        }
    }        
}