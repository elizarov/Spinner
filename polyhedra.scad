// shape height
//sh = 30;
// border width
//bw = 3;
// connector width
//cw = 1.5;

//poly_fill(truncated(icosidodecahedron));
//poly_wire(tetrahedron);
//poly_wire_dual(cube);

// comparison epsilon
c_eps = 1e-6;

function sum0(v, i, r) = i < len(v) ? sum0(v, i + 1, r + v[i]) : r;
function sum(v) = sum0(v, 1, v[0]);
function avg(v) = sum(v) / len(v);
function vnorm(v) = v / norm(v);    

function sorted(a) = len(a) == 0 ? [] : 
    let(pivot   = a[floor(len(a) / 2)])
    let(lesser  = [for (y = a) if (y  < pivot) y])
    let(equal   = [for (y = a) if (y == pivot) y])
    let(greater = [for (y = a) if (y  > pivot) y])
    concat(sorted(lesser), equal, sorted(greater));    

function approx_contains(v, a) =
    let(m = [for (x = v) if (abs(x - a) < c_eps) 1])
        m != [];
       
function distinct(v, i = 0, r = []) =
    i == len(v) ? r :
        let(a = v[i])
            approx_contains(r, a) ? 
                distinct(v, i + 1, r) :
                distinct(v, i + 1, concat(r, [a]));

function face_equation(poly, fid) =
    let(vs = poly[0])
    let(fs = poly[1])
    let(f = fs[fid])
    let(a = vs[f[0]])
    let(b = vs[f[1]])
    let(c = vs[f[2]])
    let(v1 = b - a)
    let(v2 = c - a)
    let(n = vnorm(cross(v1, v2)))
    [n.x, n.y, n.z, n * a];

function face_equations(poly) =
    let(fs = poly[1])
    [for(fid = [0:len(fs) - 1]) face_equation(poly, fid)];
    
function inradius(poly, face_equations = undef) =
    let(fe = is_undef(face_equations) ? face_equations(poly) : face_equations)
    distinct([
        for (e = fe)
            abs(e[3])
    ]);
    
function midradius(poly) = 
    let(vs = poly[0])
    let(fs = poly[1])
    distinct([
        for (f = fs)
            for (i = [0:len(f) - 1])
                let(j = (i + 1) % len(f))
                    norm((vs[f[i]] + vs[f[j]]) / 2)
    ]);
            
function circumradius(poly) = 
    let(vs = poly[0])
    distinct([
        for (v = vs) 
            norm(v)
    ]); 
     
function is_adjacent_face(f1, f2) = 
    len([for(u = f1) for(v = f2) if (u == v) 1]) == 2;
        
function dihedral_angle(poly) = 
    let(vs = poly[0])
    let(fs = poly[1])
    let(fe = face_equations(poly))
    distinct([
        for(k = [0:len(fs) - 2])
            let(f1 = fs[k])
            let(e1 = fe[k])
            for(t = [k + 1:len(fs) - 1])
                let(f2 = fs[t])
                if (is_adjacent_face(f1, f2))
                    let(e2 = fe[t])
                        180 - acos(e1.x * e2.x + e1.y * e2.y + e1.z * e2.z)
    ]);
        
function face_center(poly, fid = 0) =
    let(eq = face_equation(poly, fid))
    let(d = eq[3])
    [eq.x * d, eq.y * d, eq.z * d];
        
function face_dist(poly, fid = 0) =
    norm(face_center(poly, fid));        
    
function diameter(poly, fid = 0) =
    let(c = face_center(poly, fid))
    let(vs = poly[0])
    // check if this polyhedron has a vertex opposing a center of a face
    let(ff = [for (p = vs) if (c * p < 0 && norm(cross(c, p)) < c_eps) p])
    // no -- diameter face to face    
    ff == [] ? 2 * norm(c) :
        // yes -- diameter face to this vertex
        norm(c - ff[0]);
   
function filter_face_ids(poly, fill_reg) = 
    let(fs = poly[1])
    [
        for(k = fill_reg)
            for (i = [0:len(fs) - 1]) 
                if (len(fs[i]) == k)
                    i
    ];
    
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
    
// ------------------- Dual -------------------
 
function dual(poly, face_equations = undef) =
    let(fe = is_undef(face_equations) ? face_equations(poly) : face_equations)
    let(mrs = midradius(poly))
    let(mr = mrs[0])    
    let(vs = poly[0])
    let(fs = poly[1])
    // dual vertices coords
    let(dvs = [
        for (e = fe) 
            let(s = mr * mr / e[3]) 
                [e.x * s, e.y * s, e.z * s]
    ])
    // dual faces
    let(dfs = [
        for (i = [0:len(vs) - 1]) 
            sort_face(dvs, [
                for (j = [0:len(fs) - 1])
                    if (search(i, fs[j]) != [])
                        j
            ])
    ])
        [dvs, dfs];
                    
// ------------------- Truncation -------------------                    
                    
function find_vp(vp, u, v) = 
    search([[u, v]], vp)[0];    

// makes regular faces remain regular    
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

// ------------------- Rectification -------------------
                        
function find_vpr(vp, u, v) =
    u < v ? find_vp(vp, u, v) : find_vp(vp, v, u);

function rectified(poly) =
    let(vs = poly[0])
    let(fs = poly[1])
    // pair of original vertices
    let(vp = [
        for (f = fs) 
            for (i = [0:len(f) - 1])
                let(j = (i + 1) % len(f))
                    if (f[i] < f[j])
                        [f[i], f[j]]
    ])
    // truncated vertices coords
    let(rvs = [
        for (p = vp) 
            let(a = vs[p[0]])
            let(b = vs[p[1]])
                (a + b) / 2
    ])
    // faces from the original faces
    let(rf1 = [
        for (f = fs) [
            for (i = [0:len(f) - 1])
                let(j = (i + 1) % len(f))
                    find_vpr(vp, f[i], f[j])
        ]
    ]) 
    // faces from the original vertices
    let(rf2 = [
        for (i = [0:len(vs) - 1]) 
            sort_face(rvs, [
                for (f = fs) 
                    for (j = [0:len(f) - 1])
                        if (f[j] == i)
                            find_vpr(vp, i, f[(j + 1) % len(f)])
            ])
    ])
        [rvs, concat(rf1, rf2)];

// ------------------- Cantellation -------------------

// assuming faces were regular, ensures that resulting 4-faces are regular
function cantellataion_frac(poly) = 
    let(fs = poly[1])
    let(n = len(fs[0])) // primary face size (assuming it is regular)                    
    let(fa = 380 / n)
    let(da = dihedral_angle(poly)[0]) // angle between primary faces
        1 / (sin(da / 2) / tan(fa / 2) + 1);

function find_fv(fv, k, v) = 
    search([[k, v]], fv)[0];    

function cantellated(poly, cf = undef) =
    let(vs = poly[0])
    let(fs = poly[1])
    let(fe = face_equations(poly))
    let(c0 = is_undef(cf) ? cantellataion_frac(poly) : cf)
    // new vertex ids -- tuples of original faces and vertices
    let(fv = [
        for (k = [0:len(fs) - 1])
            for (v = fs[k])
                [k, v]
    ])
    // next vertices coords (old face, shifted towards center)
    let(cvs = [
        for (k = [0:len(fs) - 1])
            let (f = fs[k]) 
            let (eq = fe[k])
            let(d = eq[3])
            let(c = [eq.x * d, eq.y * d, eq.z * d]) // projected face center
            for (v = f) 
                let(a = vs[v])
                    a + (c - a) * c0
    ])
    // faces from the original faces
    let(cf1 = [
        for (k = [0:len(fs) - 1])
            [for (v = fs[k]) find_fv(fv, k, v)]
    ]) 
    // faces from the original vertices
    let(cf2 = [
        for (v = [0:len(vs) - 1]) 
            sort_face(cvs, [
                for (k = [0:len(fs) - 1])
                    for (vv = fs[k])
                        if (vv == v)
                            find_fv(fv, k, v)
            ])
    ])
    // faces from the original edges
    let(cf3 = [
        for (k1 = [0:len(fs) - 1]) 
            let(f1 = fs[k1])
            for (i1 = [0:len(f1) - 1])
                let(u = f1[i1])
                let(v = f1[(i1 + 1) % len(f1)])
                if (u < v) 
                    let(k1u = find_fv(fv, k1, u))
                    let(k1v = find_fv(fv, k1, v))
                    for (k2 = [0:len(fs) - 1])
                        let(f2 = fs[k2])
                        for (i2 = [0:len(f2) - 1])
                            if (f2[i2] == v && f2[(i2 + 1) % len(f2)] == u)
                                let(k2u = find_fv(fv, k2, u))
                                let(k2v = find_fv(fv, k2, v))
                                    [k1v, k1u, k2u, k2v]            
    ])
        [cvs, concat(cf1, cf2, cf3)];


    
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

module poly_fill(
    poly, sh = sh, fid = 0
) {
    validate(poly);
    pd = diameter(poly, fid);
    translate([0, 0, sh * face_dist(poly, fid) / pd])
        face_rotate(poly, fid)
            scale(sh / pd) 
                polyhedron(poly[0], poly[1]);    
}

module poly_fill0(
    poly, pd, 
    bw = bw, fid = 0
) {
    scale(bw / pd) 
        polyhedron(poly[0], poly[1]);    
}

module poly_wire_edges_impl(
    poly, pd, bw, s, 
    wire_poly, fid,
    fill_edges = [],
    eps = 0,
) {
    vs = poly[0];
    fs = poly[1];
    for (f = fs) {
        n = len(f);
        for (i = [0:n - 1]) {
            j = (i + 1) % n;
            p = f[i];
            q = f[j];
            if (p < q && search([[p, q]], fill_edges) == [[]]) {
                p0 = vs[p] * s; 
                q0 = vs[q] * s;  
                pq = vnorm(q0 - p0) * eps;
                p1 = p0 + pq;
                q1 = q0 - pq; 
                hull() {
                    translate(p1) poly_fill0(wire_poly, pd, bw, fid);
                    translate(q1) poly_fill0(wire_poly, pd, bw, fid);
                }                
            }
        }
    }   
}

module poly_wire_faces_impl(
    poly, pd, bw, s, fid, fill_all
) {
    vs = poly[0];
    fs = poly[1];
    module rec_hull(v, i) {
        if (i == 0) {
            translate(v[0]) poly_fill0(poly, pd, bw, fid);
        } else {
            hull() {
                translate(v[i]) poly_fill0(poly, pd, bw, fid);
                rec_hull(v, i - 1);
            }
        }
    }
    for (k = fill_all) {
        v = [for (p = fs[k]) vs[p] * s];
        rec_hull(v, len(v) - 1);    
    }  
}

function fill_all(poly, fill, fill_reg) =
    concat(fill, filter_face_ids(poly, fill_reg)); 

function fill_edges(poly, fill_all) =
    let(fs = poly[1])
    [
        for (k = fill_all) 
            let(f = fs[k])
            for (i = [0:len(f) - 1])
                let(p = f[i])
                let(q = f[(i + 1) % len(f)])
                    p < q ? [p, q] : [q, p]
    ];
            
function fill_verts(poly, fill_all) =
    let(fs = poly[1])
    [
        for (k = fill_all) 
            let(f = fs[k])
            for (p = f)
                p
    ];        

// fill - a list of face ids to fill
// fill_reg - a list of face sizes to fill
module poly_wire0(
    poly, pd,
    sh = sh, bw = bw, fid = 0, 
    fill = [], fill_reg = [],
    eps = 0
) {
    s = (sh - bw) / pd;
    vs = poly[0];
    fs = poly[1];
    fill_all = fill_all(poly, fill, fill_reg);
    fill_edges = fill_edges(poly, fill_all);
    fill_verts = fill_verts(poly, fill_all); 
    // filled faces
    poly_wire_faces_impl(poly, pd, bw, s, fid, fill_all); 
    // edges
    poly_wire_edges_impl(poly, pd, bw, s, poly, fid, fill_edges, eps);
    // vertices
    if (eps != 0) {
        for (i = [0:len(vs) - 1]) {
            if (search(i, fill_verts) == []) {
                p0 = vs[i] * s;
                translate(p0) poly_fill0(poly, pd, bw, fid); 
            }
        }
    }
}

module poly_wire(
    poly, 
    sh = sh, bw = bw, fid = 0, 
    fill = [], fill_reg = [], 
    eps = 0
) {
    validate(poly);
    pd = diameter(poly, fid);
    translate([0, 0, sh * face_dist(poly, fid) / pd])
        face_rotate(poly, fid)
            poly_wire0(poly, pd, sh, bw, fid, fill, fill_reg, eps);
}

module poly_wire_dual0(
    poly, pd, 
    sh = sh, bw = bw, cw = cw, dw = undef, fid = 0,
    fill = [], fill_reg = [],
) {
    dual = dual(poly);
    dradius = circumradius(dual)[0];

    dw0 = is_undef(dw) ? bw : dw;
    
    s = (sh - bw) / pd;
    dscale = (sh - dw0) / (2 * dradius);
    dfactor = dscale / s;
    
    vs = poly[0];
    fs = poly[1];
    dvs = dual[0];
    fill_all = fill_all(poly, fill, fill_reg);
    fill_edges = fill_edges(poly, fill_all);
    poly_wire_faces_impl(poly, pd, bw, s, fid, fill_all); 
    poly_wire_edges_impl(poly, pd, bw, s, poly, fid, fill_edges = fill_edges);
    poly_wire_edges_impl(dual, pd, dw0, dscale, poly, fid);
    // wires between primary and dual
    if (cw != 0) {
        cs = (sh - cw) / pd;
        for (i = [0:len(fs) - 1]) {
            c0 = dvs[i] * cs * dfactor;
            for (u = fs[i]) {
                p0 = vs[u] * cs;
                hull() {
                    translate(c0) poly_fill0(poly, pd, cw, fid);
                    translate(p0) poly_fill0(poly, pd, cw, fid);
                }       
            }
        }
    }
}

module poly_wire_dual(
    poly, 
    sh = sh, bw = bw, cw = cw, dw = undef, fid = 0,
    fill = [], fill_reg = [],
) {
    validate(poly);
    pd = diameter(poly, fid);
    translate([0, 0, sh * face_dist(poly, fid) / pd])
        face_rotate(poly, fid) 
            poly_wire_dual0(poly, pd, sh, bw, cw, dw, fid, fill, fill_reg);
}

module validate(poly) {
    vs = poly[0];
    fs = poly[1];
    face_sizes = distinct([
        for (f = fs) len(f)
    ]);
    echo("# Vertices:", N = len(vs));
    echo("# Faces:", N = len(fs), face_sizes = face_sizes);
    if (len(face_sizes) > 1) {
        for (m = face_sizes) {
            echo(m, "-faces:", N = len([for (f = fs) if (len(f) == m) 1]));
        }    
    }
    fe = face_equations(poly);
    inradius = inradius(poly, fe);
    midradius = midradius(poly);
    circumradius = circumradius(poly);
    dihedral_angle = dihedral_angle(poly);
    echo(inradius = inradius);
    echo(midradius = midradius);
    echo(circumradius = circumradius);
    echo(dihedral_angle = dihedral_angle);
    for (i = [0:len(fs) - 1]) {
        f = fs[i];
        e = fe[i];
        c = avg([for (u = f) vs[u]]); // face center
        for (j = [0:len(f) - 1]) {
            p = vs[f[j]];
            dist = abs(concat(p, [-1]) * e);
            assert(dist < c_eps, ["Face is not planar #", i, f, c]);
            q = vs[f[(j + 1) % len(f)]];
            dir = concat(cross(vnorm(q - c), vnorm(p - c)), [0]) * e;
            assert(dir < c_eps, ["Face is not clockwise #", i, f, c]);
        }
    }
}

// --------------------- 5 Platonic Solids ---------------------

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

// --------------------- 13 Arhimedean Solids ---------------------

truncated_tetrahedron = truncated(tetrahedron);
cuboctahedron = rectified(cube);
truncated_cube = truncated(cube);
truncated_octahedron = truncated(octahedron);
rhombicuboctahedron = cantellated(cube);
// rhombitruncated_cuboctahedron
// snub_cube
icosidodecahedron = rectified(dodecahedron);
truncated_dodecahedron = truncated(dodecahedron);
truncated_icosahedron = truncated(icosahedron);
rhombicosidodecahedron = cantellated(dodecahedron);
// rhombitruncated_icosidodecahedron
// snub_dodecahedron

// --------------------- 13 Catalan Solids (Arhimedean Duals) ---------------------

triakis_tetrahedron = dual(truncated_tetrahedron);
rhombic_dodecahedron = dual(cuboctahedron);
triakis_octahedron = dual(truncated_cube);
tetrakis_hexahedron = dual(truncated_octahedron);
deltoidal_icositetrahedron = dual(rhombicuboctahedron);
// disdyakis_dodecahedron = dual(rhombitruncated_cuboctahedron);
// pentagonal_icositetrahedron = dual(snub_cube);
rhombic_triacontahedron = dual(icosidodecahedron);
triakis_icosahedron = dual(truncated_dodecahedron);
pentakis_dodecahedron = dual(truncated_icosahedron);
deltoidal_hexecontahedron = dual(rhombicosidodecahedron);
// disdyakis_triacontahedron = dual(rhombitruncated_icosidodecahedron);
// pentagonal_hexecontahedron = dual(snub_dodecahedron);
