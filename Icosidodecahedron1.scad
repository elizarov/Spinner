// shape height
sh = 30;
// border width
bw = 3;

translate([0, 0, sh / 2])
    idface_rotate()
         icosidodecahedron1();

// ------------ dodecahedron ------------

phi = (sqrt(5) + 1) / 2;

// dodecahedron vertices

dvert = [
    [-1/phi, 0, -phi], // 0
    [1/phi, 0, -phi], // 1
    [-1, -1, -1], // 2
    [-1, 1, -1], // 3
    [1, -1, -1], // 4
    [1, 1, -1], // 5
    [0, -phi, -1/phi], // 6
    [0, phi, -1/phi], // 7
    [-phi, -1/phi, 0], // 8
    [-phi, 1/phi, 0], // 9
    [phi, -1/phi, 0], // 10
    [phi, 1/phi, 0], // 11
    [0, -phi, 1/phi], // 12
    [0, phi, 1/phi], // 13
    [-1, -1, 1], // 14
    [-1, 1, 1], // 15
    [1, -1, 1], // 16
    [1, 1, 1], // 17
    [-1/phi, 0, phi], // 18
    [1/phi, 0, phi] // 19    
];

dface = [
    [0, 1, 5, 7, 3],
    [1, 0, 2, 6, 4],
    [0, 3, 9, 8, 2],
    [1, 4, 10, 11, 5],
    [8, 14, 12, 6, 2],
    [6, 12, 16, 10, 4],
    [7, 13, 15, 9, 3],
    [13, 7, 5, 11, 17],
    [8, 9, 15, 18, 14],
    [11, 10, 16, 19, 17],
    [18, 19, 16, 12, 14],
    [19, 18, 15, 13, 17]
];

// icosidodecahedron vertices -- rectified dodecahedron

idrectify = [
    for (f = dface) 
        for (i = [0:4]) 
            let(j = (i + 1) % 5)
            let(p = f[i])
            let(q = f[j])
            if (p < q) 
                [p, q]
];

idvert = [
    for (r = idrectify)
        (dvert[r[0]] + dvert[r[1]]) / 2
];
       
function idvn(p, q) =
    p > q ? idvn(q, p) :
        search([[p, q]], idrectify)[0]; 
        
// pentagonal faces    
idface_p = [
    for (f = dface) [
        for (i = [0:4]) 
            let(j = (i + 1) % 5)
                idvn(f[i], f[j])
    ]
];            
            
// triangual faces
idface_t = [
    for (i = [0:19]) [
        for (f = dface) 
            for(j = [0:4]) 
                if (f[j] == i) 
                    idvn(i, f[(j + 1) % 5])
    ]
];

idface = concat(idface_p, idface_t);       

idf0 = idface_p[0];
idf0c = (idvert[idf0[0]] + idvert[idf0[1]] + idvert[idf0[2]] + idvert[idf0[3]] + idvert[idf0[4]]) / 5;
idfr = norm(idf0c);
idfa = atan(idf0c.y / idf0c.z);
                
module icosidodecahedron1() {
    s = (sh - bw) / 2 / idfr;
    for (f = idface_p) {
        for (i = [0:4]) {
            j = (i + 1) % 5;
            p = f[i];
            q = f[j];
            hull() {
                translate(idvert[p] * s) icosidodecahedron0();
                translate(idvert[q] * s) icosidodecahedron0();
            }                
        }
    }   
}

module icosidodecahedron0() {
    scale(bw / 2 / idfr)    
        polyhedron(idvert, idface);
}

module idface_rotate() {
    rotate([idfa, 0, 0]) children();
}