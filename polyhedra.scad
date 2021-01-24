// shape height
//sh = 30;
// border width
//bw = 3;

//poly_wire(octahedron);

function sum0(v, i, r) = i < len(v) ? sum0(v, i + 1, r + v[i]) : r;
function sum(v) = sum0(v, 1, v[0]);
function avg(v) = sum(v) / len(v);

function face_center(poly, fid = 0) =
    let(vs = poly[0])
    let(f = poly[1][fid])
    avg([for (p = f) vs[p]]);
    
function face_dist(poly, fid = 0) =
    norm(face_center(poly, fid));        
    
module face_rotate(poly, fid = 0) {
    c = face_center(poly, fid);
    az = c.x == 0 ? 90 : -atan(c.y / c.x);
    rz = [
        [cos(az), -sin(az), 0],
        [sin(az), cos(az), 0],
        [0, 0, 1],
    ] * c;
    ay = atan(rz.x / rz.z);
    rotate([0, ay, 0]) rotate([0, 0, az])
        children();
}

module poly_fill0(poly, bw = bw, fid = 0) {
    scale(bw / 2 / face_dist(poly, fid)) 
        polyhedron(poly[0], poly[1]);    
}

module poly_wire0(poly, sh = sh, bw = bw, fid = 0) {
    s = (sh - bw) / 2 / face_dist(poly, fid);
    vs = poly[0];
    for (f = poly[1]) {
        n = len(f);
        for (i = [0:n - 1]) {
            j = (i + 1) % n;
            p = f[i];
            q = f[j];
            if (p < q) {
                hull() {
                    translate(vs[p] * s) poly_fill0(poly, bw, fid);
                    translate(vs[q] * s) poly_fill0(poly, bw, fid);
                }                
            }
        }
    }        
}

module poly_wire(poly, sh = sh, bw = bw, fid = 0) {
    translate([0, 0, sh / 2])
        face_rotate(poly, fid)
            poly_wire0(poly, sh, bw, fid);
}

// --------------------- Platonic Solids ---------------------

phi = (sqrt(5) + 1) / 2;

tetrahedron = [[
    [-1, 0, -1/sqrt(2)], // 0
    [1, 0, -1/sqrt(2)], // 1
    [0, -1, 1/sqrt(2)], // 2
    [0, 1, 1/sqrt(2)] // 3
], [
    [0, 1, 3],
    [0, 2, 1],
    [0, 3, 2],
    [1, 2, 3]
]];

cube = [[
    [1, 1, -1], // 0
    [-1, 1, -1], // 1
    [-1, -1, -1], // 2
    [1, -1, -1], // 3
    [1, 1, 1], // 4
    [-1, 1, 1], // 5
    [-1, -1, 1], // 6
    [1, -1, 1], // 7
], [
    [0, 1, 2, 3],
    [0, 4, 5, 1],
    [1, 5, 6, 2],
    [2, 6, 7, 3],
    [3, 7, 4, 0],
    [4, 7, 6, 5]
]];

octahedron = [[
    [0, 0, -1], // 0
    [1, 0, 0], // 1
    [0, 1, 0], // 2
    [-1, 0, 0], // 3
    [0, -1, 0], // 4
    [0, 0, 1] // 5
], [
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 4],
    [0, 4, 1],
    [5, 2, 1],
    [5, 3, 2],
    [5, 4, 3],
    [5, 1, 4]
]];

dodecahedron = [[
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
], [
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
]];

icosahedron = [[
    [0, -1, -phi], // 0
    [0, 1, -phi], // 1
    [-phi, 0, -1], // 2
    [phi, 0, -1], // 3
    [-1, -phi, 0], // 4
    [-1, phi, 0], // 5
    [1, -phi, 0], // 6
    [1, phi, 0], // 7
    [-phi, 0, 1], // 8
    [phi, 0, 1], // 9
    [0, -1, phi], // 10
    [0, 1, phi] // 11    
], [
    [0, 1, 2],
    [1, 0, 3],
    [0, 2, 4],
    [2, 1, 5],
    [1, 3, 7],
    [3, 0, 6],
    [1, 7, 5],
    [0, 4, 6],
    [2, 8, 4],
    [2, 5, 8],
    [3, 6, 9],
    [3, 9, 7],
    [4, 10, 6],
    [5, 7, 11],
    [8, 10, 4],
    [5, 11, 8],
    [9, 11, 7],
    [6, 10, 9],
    [8, 11, 10],
    [9, 10, 11]
]];
