// shape height
//sh = 30;
// border width
//bw = 3;

//poly_wire(truncated(cube), fill = [0]);

function sum0(v, i, r) = i < len(v) ? sum0(v, i + 1, r + v[i]) : r;
function sum(v) = sum0(v, 1, v[0]);
function avg(v) = sum(v) / len(v);

function face_center(poly, fid = 0) =
    let(vs = poly[0])
    let(f = poly[1][fid])
    avg([for (p = f) vs[p]]);
    
function face_dist(poly, fid = 0) =
    norm(face_center(poly, fid));        
    
function find_vp(vp, u, v) = 
    search([[u, v]], vp)[0];    
    
function vnorm(v) = v / norm(v);    
          
function next_sorted_face(vs, f, c, u) =
    let(a = vnorm(vs[u] - c))
    let(m = [
        for (v = f) 
            v == u ? -3 :
            let(b = vnorm(vs[v] - c))
                cross(a, b) * c >= 0 ? -2 :
                    a * b
    ]) 
    let(mm = max(m))
        f[search(mm, m)[0]];
        
function sort_face0(vs, f, c, i, r) =
    i == len(f) - 1 ? r : 
        sort_face0(vs, f, c, i + 1, 
            concat(r, [next_sorted_face(vs, f, c, r[i])])
        );
    
function sort_face(vs, f) =
    let(c = avg([for (k = f) vs[k]]))
        sort_face0(vs, f, c, 0, [f[0]]);
    
function truncation_frac(n) =
    let(beta = 90 - 180/n)
        1 / (2 * (1 + sin(beta)));
     
function truncated(poly, tr = undef) =
    let(vs = poly[0])
    let(fs = poly[1])
    let(t0 = is_undef(tr) ? truncation_frac(len(fs[0])) : tr)
    // pair of original vertices
    let(vp = [
        for (f = fs) 
            for (i = [0:len(f) - 1])
                let(j = (i + 1) % len(f))
                    [f[i], f[j]]
    ])
    // truncated vertices coords
    let(tvs = [
        for (p = vp) 
            let(a = vs[p[0]])
            let(b = vs[p[1]])
                a + t0 * (b - a)
    ])
    // faces from the original faces
    let(tf1 = [
        for (f = fs) [
            for (t = [0:2 * len(f) - 1])
                let(i = t / 2)
                let(j = (i + 1) % len(f))
                    t % 2 == 0 ? find_vp(vp, f[i], f[j]) : find_vp(vp, f[j], f[i])
        ]
    ]) 
    // faces from the original vertices
    let(tf2 = [
        for (i = [0:len(vs) - 1]) 
            sort_face(tvs, [
                for (f = fs) 
                    for (j = [0:len(f) - 1])
                        if (f[j] == i)
                            find_vp(vp, i, f[(j + 1) % len(f)])
            ])
    ])
        [tvs, concat(tf1, tf2)];
    
module face_rotate(poly, fid = 0) {
    c = face_center(poly, fid);
    az = c.x == 0 ? 90 : -atan(c.y / c.x);
    rz = [
        [cos(az), -sin(az), 0],
        [sin(az), cos(az), 0],
        [0, 0, 1],
    ] * c;
    ay = -atan(rz.x / rz.z);
    rotate([0, ay, 0]) rotate([0, 0, az])
        children();
}

module poly_fill0(poly, bw = bw, fid = 0) {
    scale(bw / 2 / face_dist(poly, fid)) 
        polyhedron(poly[0], poly[1]);    
}

module poly_wire0(poly, sh = sh, bw = bw, fid = 0, fill = []) {
    s = (sh - bw) / 2 / face_dist(poly, fid);
    vs = poly[0];
    fs = poly[1];
    module rec_hull(v, i) {
        if (i == 0) {
            translate(v[0]) poly_fill0(poly, bw, fid);
        } else {
            hull() {
                translate(v[i]) poly_fill0(poly, bw, fid);
                rec_hull(v, i - 1);
            }
        }
    }
    for (k = fill) {
        v = [for (p = fs[k]) vs[p] * s];
        rec_hull(v, len(v) - 1);    
    }
    fill_edges = [
        for (k = fill) 
            let(f = fs[k])
            for (i = [0:len(f) - 1])
                let(p = f[i])
                let(q = f[(i + 1) % len(f)])
                    p < q ? [p, q] : [q, p]
    ];
    for (f = fs) {
        n = len(f);
        for (i = [0:n - 1]) {
            j = (i + 1) % n;
            p = f[i];
            q = f[j];
            if (p < q && search([[p, q]], fill_edges) == [[]]) {
                hull() {
                    translate(vs[p] * s) poly_fill0(poly, bw, fid);
                    translate(vs[q] * s) poly_fill0(poly, bw, fid);
                }                
            }
        }
    }        
}

module poly_wire(poly, sh = sh, bw = bw, fid = 0, fill = []) {
    translate([0, 0, sh / 2])
        face_rotate(poly, fid)
            poly_wire0(poly, sh, bw, fid, fill);
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

// --------------------- Arhimedian Solids ---------------------

truncated_octahedron = truncated(octahedron);

