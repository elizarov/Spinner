// shape height
//sh = 30;
// border width
//bw = 3;
// connector width
//cw = 1.5;

//poly_fill(great_stellated_dodecahedron);
//poly_fill(icosahedron);
//poly_wire(tetrahedron);
//poly_wire_dual(cube);

// comparison epsilon
c_eps = 1e-10;

// snub computation precision
s_eps = 1e-15;

function sum0(v, i, r) = i < len(v) ? sum0(v, i + 1, r + v[i]) : r;
function sum(v) = sum0(v, 1, v[0]);
function avg(v) = sum(v) / len(v);
function vnorm(v) = 
    let(l = norm(v))
    l < c_eps ? v : v / l;    

function sqr(x) = x * x;

function sorted(a) = len(a) == 0 ? [] : 
    let(pivot   = a[floor(len(a) / 2)])
    let(lesser  = [for (y = a) if (y  < pivot) y])
    let(equal   = [for (y = a) if (y == pivot) y])
    let(greater = [for (y = a) if (y  > pivot) y])
    concat(sorted(lesser), equal, sorted(greater));    

function approx_contains(v, a) =
    let(m = [for (x = v) if (abs(x - a) < c_eps) 1])
        len(m) > 0;
       
function distinct(v, i = 0, r = []) =
    i == len(v) ? r :
        let(a = v[i])
            approx_contains(r, a) ? 
                distinct(v, i + 1, r) :
                distinct(v, i + 1, concat(r, [a]));
    
// plane equation by 3 point    
function plane3(a, b, c) =    
    let(v1 = b - a)
    let(v2 = c - a)
    let(n = vnorm(cross(v1, v2)))
    [n.x, n.y, n.z, n * a];
    
function is_outside_plane(poly, eq) = 
    let(n = [eq.x, eq.y, eq.z])
    norm(n) < c_eps ? false :
        let(vs = poly[0])  
        let(ds = [for(a = vs) n * a - eq[3]])
        let(pos = [for(d = ds) if (d > c_eps) 1])
        let(neg = [for(d = ds) if (d < -c_eps) 1])
            len(pos) == 0 || len(neg) == 0;
    
function find_outside_plane(poly) =
    let(vs = poly[0])    
    let(ops = [
        for(i = [0:len(vs) - 3])
            for(j = [i + 1:len(vs) - 2])
                for(k = [j + 1:len(vs) - 1])
                    let(eq = plane3(vs[i], vs[j], vs[k]))
                    if (is_outside_plane(poly, eq))
                        eq
    ])
    assert(len(ops) > 0, "Cannot find outside plane")
        ops[0];
                            
function default_face_equation(poly) =
    let(eq0 = face_equation(poly, 0))
    is_outside_plane(poly, eq0) ? eq0 :
        find_outside_plane(poly);    

function face_equation(poly, fid = undef) =
    is_undef(fid) ? default_face_equation(poly) :
    let(vs = poly[0])
    let(fs = poly[1])
    let(f = fs[fid])
    plane3(vs[f[0]], vs[f[1]], vs[f[2]]);

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
        
function face_center(poly, fid = undef) =
    let(eq = face_equation(poly, fid))
    let(d = eq[3])
    [eq.x * d, eq.y * d, eq.z * d];
        
function face_dist(poly, fid = undef) =
    norm(face_center(poly, fid));        
    
function diameter(poly, fid = undef) =
    let(c = face_center(poly, fid))
    let(vs = poly[0])
    // check if this polyhedron has a vertex opposing a center of a face
    let(ff = [for (p = vs) if (c * p < 0 && norm(cross(c, p)) < c_eps) p])
    // no -- diameter face to face    
    ff == [] ? 2 * norm(c) :
        // yes -- diameter face to this vertex
        norm(c - ff[0]);
   
function filter_face_ids(poly, fill_face_sizes) = 
    let(fs = poly[1])
    [
        for(k = fill_face_sizes)
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
        1 / (1 + sin(beta));
     
// tf -- truncation fraction from 0 to 1                    
//       0 - no truncation
//       1 - full truncation (rectification)                    
function truncated(poly, tf = undef) =
    tf <= 0 ? poly : tf >= 1 ? rectified(poly) :                    
    let(vs = poly[0])
    let(fs = poly[1])
    let(t0 = is_undef(tf) ? truncation_frac(len(fs[0])) : tf)
    // pairs of original vertices
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
                a + t0 * (b - a) / 2
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
    let(fa = 360 / n)
    let(da = dihedral_angle(poly)[0]) // angle between primary faces
        1 / (sin(da / 2) / tan(fa / 2) + 1);

function find_fv(fv, k, v) = 
    search([[k, v]], fv)[0];    

// cf -- cantellation fraction from 0 to 1                    
//       0 - no cantellation
//       1 - full cantellation (dual)       
function cantellated(poly, cf = undef) =
    cf <= 0 ? poly : cf >= 1 ? dual(poly) :
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

// ------------------- Bevel -------------------

function find_vpf(vpf, u, v, k) = 
    search([[u, v, k]], vpf)[0];    

function bevel_frac(poly) = 
    let(fs = poly[1])
    let(n = len(fs[0])) // primary face size (assuming it is regular)                    
    let(fa = 360 / n)
    let(da = dihedral_angle(poly)[0]) // angle between primary faces
    let(tf = truncation_frac(n))
    let(cf = (1 - tf) / (sin(da / 2) / tan(fa / 2) + 1 - tf))
        [cf, tf];

// bf -- pair [cf, tf]
//       cf - cantellation fraction from 0 to 1                    
//       tf - truncation fraction from 0 to 1
function beveled(poly, bf = undef) =
    let(vs = poly[0])
    let(fs = poly[1])
    let(fe = face_equations(poly))
    let(b0 = is_undef(bf) ? bevel_frac(poly) : bf)
    let(c0 = b0[0])
    let(t0 = b0[1])
    // triples - pairs of original vertices both ways, face ids
    let(vpf = [
        for (k = [0:len(fs) - 1])
            let (f = fs[k]) 
            for (t = [0:2 * len(f) - 1])
                let(i = t / 2)
                let(j = (i + 1) % len(f))
                let(u = t % 2 == 0 ? f[i] : f[j])
                let(v = t % 2 == 0 ? f[j] : f[i])
                     [u, v, k]
    ])
    // beveled vertices coords
    let(bvs = [
        for (k = [0:len(fs) - 1])
            let (f = fs[k]) 
            let (eq = fe[k])
            let(d = eq[3])
            let(c = [eq.x * d, eq.y * d, eq.z * d]) // projected face center
            for (t = [0:2 * len(f) - 1])
                let(i = t / 2)
                let(j = (i + 1) % len(f))
                let(u = t % 2 == 0 ? f[i] : f[j])
                let(v = t % 2 == 0 ? f[j] : f[i])
                let(a = vs[u])
                let(b = vs[v])
                let(ac = a + (c - a) * c0) // cantellated a (shifted towards center)
                let(bc = b + (c - b) * c0) // cantellated b (shifted towards center)
                    ac + t0 * (bc - ac) / 2 // [a,b] truncation 
    ])
    // faces from the original faces (like in truncation)
    let(bf1 = [
        for (k = [0:len(fs) - 1])
            let (f = fs[k]) 
            [
                for (t = [0:2 * len(f) - 1])
                    let(i = t / 2)
                    let(j = (i + 1) % len(f))
                        t % 2 == 0 ? find_vpf(vpf, f[i], f[j], k) : find_vpf(vpf, f[j], f[i], k)
            ]
    ]) 
    // faces from the original vertices
    let(bf2 = [
        for (v = [0:len(vs) - 1]) 
            sort_face(bvs, [
                for (k = [0:len(fs) - 1])
                    let(f = fs[k])
                    for (t = [0:2 * len(f) - 1])
                        let(i = t / 2)
                        if (f[i] == v)
                            let(j = (i + len(f) + (t % 2 == 0 ? -1 : 1)) % len(f))
                                find_vpf(vpf, v, f[j], k)
            ])
    ])
    // faces from the original edges (like in cantellation)
    let(bf3 = [
        for (k1 = [0:len(fs) - 1]) 
            let(f1 = fs[k1])
            for (i1 = [0:len(f1) - 1])
                let(u = f1[i1])
                let(v = f1[(i1 + 1) % len(f1)])
                if (u < v) 
                    let(k1u = find_vpf(vpf, u, v, k1))
                    let(k1v = find_vpf(vpf, v, u, k1))
                    for (k2 = [0:len(fs) - 1])
                        let(f2 = fs[k2])
                        for (i2 = [0:len(f2) - 1])
                            if (f2[i2] == v && f2[(i2 + 1) % len(f2)] == u)
                                let(k2u = find_vpf(vpf, u, v, k2))
                                let(k2v = find_vpf(vpf, v, u, k2))
                                    [k1v, k1u, k2u, k2v]            
    ])
        [bvs, concat(bf1, bf2, bf3)];

// ------------------- Snub -------------------

// a * x^2 + b * x + c = 0
function solve3(a, b, c) = 
    (-b + sqrt(sqr(b) - 4 * a * c)) / (2 * a);

// solves for snub_compute_a(da, fa, cf, sa) == rf
function snub_compute_sa(da, fa, cf) =
    let(rf = 1 - cf)
    let(cm = (1 - cos(da)) / 2)
    let(cp = (1 + cos(da)) / 2)
    let(cos_ga = solve3(cp * sqr(rf), 2 * cm * rf * cos(fa/2), -sqr(cos(fa/2)) * (cm + sqr(rf))))
        fa/2 - acos(cos_ga);

function snub_compute_a(da, fa, cf, sa) =
    let(h = 1 / (2 * tan(fa / 2)))
    let(ga = fa/2 - sa)
    let(rf = 1 - cf)
    let(t = rf / (2 * sin(fa/2)))
    norm([
        2 * t * sin(ga),
        (h - t * cos(ga)) * (cos(da) - 1),
        (h - t * cos(ga)) * sin(da)
    ]);
    
function snub_compute_b(da, fa, cf, sa) =
    let(h = 1 / (2 * tan(fa/2)))
    let(ga = fa/2 - sa)
    let(ha = fa/2 + sa)
    let(rf = 1 - cf)
    let(t = rf / (2 * sin(fa/2)))
    norm([
        t * (sin(ga) - sin(ha)), 
        (h - t * cos(ga)) * cos(da) - (h - t * cos(ha)),
        (h - t * cos(ga)) * sin(da)
    ]);

function snub_compute_cf(da, fa, cf_l = 0, cf_r = 1) = 
    let(cf = (cf_l + cf_r) / 2)
    cf_r - cf_l < s_eps ? cf :
    let(sa = snub_compute_sa(da, fa, cf))
    let(rf = 1 - cf)
    // error goes from postive to negative to NaN as cf goes from 0 to 1                            
    let(err = snub_compute_b(da, fa, cf, sa) - rf)
        err <= 0 ? 
            snub_compute_cf(da, fa, cf, cf_r) : 
            snub_compute_cf(da, fa, cf_l, cf);    

// assuming faces were regular, ensures that resulting 3-faces are regular
function snub_frac(poly) = 
    let(fs = poly[1])
    let(n = len(fs[0])) // primary face size (assuming it is regular)                    
    let(fa = 360 / n)
    let(da = dihedral_angle(poly)[0]) // angle between primary faces
    let(cf = snub_compute_cf(da, fa))
    let(sa = snub_compute_sa(da, fa, cf))
        [cf, sa];

id_mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

function cross_mat(u) = [
    [0, -u.z, u.y],
    [u.z, 0, -u.x],
    [-u.y, u.x, 0]
];
    
function outer_prod(u, v) =
    [[u.x], [u.y], [u.z]] * [v];

function rotation_matrix(u, a) =
    cos(a) * id_mat + sin(a) * cross_mat(u) + (1 - cos(a)) * outer_prod(u, u);

// sf -- [cf, sa]
//       cf - cantellation fraction
//       sa - snub angle
function snub(poly, sf = undef) =
    let(vs = poly[0])
    let(fs = poly[1])
    let(fe = face_equations(poly))
    let(s0 = is_undef(sf) ? snub_frac(poly) : sf)
    let(c0 = s0[0])
    let(sa = s0[1])
/* BEGIN debug        
    let(n = len(fs[0])) // primary face size (assuming it is regular)                    
    let(fa = 360 / n)
    let(da = dihedral_angle(poly)[0]) // angle between primary faces
    let(rf = 1 - c0)
    echo(
        "snub: ", cf=c0, sa=sa,
            err_a = snub_compute_a(da, fa, c0, sa) - rf,
            err_b = snub_compute_b(da, fa, c0, sa) - rf
    )
// END debug */
    sa == 0 ? cantellated(poly, c0) : 
    // new vertex ids -- tuples of original faces and vertices
    let(fv = [
        for (k = [0:len(fs) - 1])
            for (v = fs[k])
                [k, v]
    ])
    // next vertices coords (old face, shifted towards center, rotated)
    let(svs = [
        for (k = [0:len(fs) - 1])
            let (f = fs[k]) 
            let (eq = fe[k])
            let(d = eq[3])
            let(c = [eq.x * d, eq.y * d, eq.z * d]) // projected face center
            let(m = (1 - c0) * rotation_matrix([eq.x, eq.y, eq.z], -sa)) // rotate clockwise!
            for (v = f) 
                let(a = vs[v])
//                    echo("face #", k, a = a, "dist", norm(a - c), "dist_new", norm((a - c) * m))
                    c + (a - c) * m
    ])
    // faces from the original faces
    let(sf1 = [
        for (k = [0:len(fs) - 1])
            [for (v = fs[k]) find_fv(fv, k, v)]
    ]) 
    // faces from the original vertices
    let(sf2 = [
        for (v = [0:len(vs) - 1]) 
            sort_face(svs, [
                for (k = [0:len(fs) - 1])
                    for (vv = fs[k])
                        if (vv == v)
                            find_fv(fv, k, v)
            ])
    ])
    // faces from the original edges
    let(sf3 = [
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
                                    each [ [k1v, k1u, k2v], [k1u, k2u, k2v] ]
    ])
        [svs, concat(sf1, sf2, sf3)];


// ------------------- Augmentation -------------------


// h -- augmentation height as a multiple of edge length             
//       0 - no augmentation
function augmented(poly, h) =
    h <= 0 ? poly :
    let(vs = poly[0])
    let(fs = poly[1])
    let(fe = face_equations(poly))
    // new vertices coords: old vertices + augmented faces
    let(avs = concat(vs, [
        for(k = [0:len(fs) - 1])
            let(f = fs[k]) 
            let(eq = fe[k])
            let(n = [eq.x, eq.y, eq.z]) // face normal
            let(d = eq[3]) // face distance
            let(c = n * d) // projected face center
            let(el = norm(vs[f[0]] - vs[f[1]])) // edge length
                c - n * el * h // augment normal to face
    ]))
    // new faces from the autmented original faces
    let(afs = [
        for (k = [0:len(fs) - 1])
            let(f = fs[k])
            let(w = len(vs) + k) // autmented vertex id
            for(i = [0:len(f) - 1]) 
                let(u = f[i])
                let(v = f[(i + 1) % len(f)])
                    [u, v, w]
    ]) 
        [avs, afs];

// ------------------- Drawing -------------------
    
module face_rotate(poly, fid = undef) {
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

module poly_place(poly, sh = sh, fid = undef, place = true, circ = false) {
    if (place) {
        s = circ ? 
            sh / (2 * max(circumradius(poly))) :
            sh / diameter(poly, fid);
        translate([0, 0, s * face_dist(poly, fid)])
            face_rotate(poly, fid)
                scale(s) 
                    children();
    } else {
        children();
    }
}

// place = true to rotate the poly onto face and scale it
// circ = true to define circumsphere diameter with sh (otherwise defines resulting height)
module poly_fill(
    poly, sh = sh, fid = undef, 
    place = true, circ = false
) {
    validate(poly);
    poly_place(poly, sh, fid, place, circ) {
        polyhedron(poly[0], poly[1]);    
    }
}

module poly_fill0(poly, s, fid = undef) {
    scale(s) 
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
                    translate(p1) poly_fill0(wire_poly, bw / pd, fid);
                    translate(q1) poly_fill0(wire_poly, bw / pd, fid);
                }                
            }
        }
    }   
}

module poly_wire_faces_impl(
    poly, pd, bw, s, 
    wire_poly, fid, 
    fill_all = []
) {
    vs = poly[0];
    fs = poly[1];
    module rec_hull(v, i) {
        if (i == 0) {
            translate(v[0]) poly_fill0(wire_poly, bw / pd, fid);
        } else {
            hull() {
                translate(v[i]) poly_fill0(wire_poly, bw / pd, fid);
                rec_hull(v, i - 1);
            }
        }
    }
    for (k = fill_all) {
        v = [for (p = fs[k]) vs[p] * s];
        rec_hull(v, len(v) - 1);    
    }  
}

function fill_all(poly, fill_face_ids, fill_face_sizes) =
    concat(fill_face_ids, filter_face_ids(poly, fill_face_sizes)); 

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

// fill_face_ids - a list of face ids to fill_face_ids
// fill_face_sizes - a list of face sizes to fill_face_ids
module poly_wire0(
    poly, pd,
    sh = sh, bw = bw, fid = undef, 
    fill_face_ids = [], fill_face_sizes = [],
    eps = 0
) {
    s = (sh - bw) / pd;
    vs = poly[0];
    fs = poly[1];
    fill_all = fill_all(poly, fill_face_ids, fill_face_sizes);
    fill_edges = fill_edges(poly, fill_all);
    fill_verts = fill_verts(poly, fill_all); 
    // filled faces
    poly_wire_faces_impl(poly, pd, bw, s, poly, fid, fill_all); 
    // edges
    poly_wire_edges_impl(poly, pd, bw, s, poly, fid, fill_edges, eps);
    // vertices
    if (eps != 0) {
        for (i = [0:len(vs) - 1]) {
            if (search(i, fill_verts) == []) {
                p0 = vs[i] * s;
                translate(p0) poly_fill0(poly, bw / pd, fid); 
            }
        }
    }
}

module poly_wire(
    poly, 
    sh = sh, bw = bw, fid = undef, 
    fill_face_ids = [], fill_face_sizes = [], 
    eps = 0,
    circ = false
) {
    validate(poly);
    pd = circ ? 
        2 * max(circumradius(poly)) :
        diameter(poly, fid);
    translate([0, 0, sh * face_dist(poly, fid) / pd])
        face_rotate(poly, fid)
            poly_wire0(poly, pd, sh, bw, fid, fill_face_ids, fill_face_sizes, eps);
}

module poly_wire_dual0(
    poly, pd, 
    sh = sh, bw = bw, cw = cw, dw = undef, fid = undef,
    fill_face_ids = [], fill_face_sizes = [], dual_fill_face_sizes = []
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
    // primary
    fill_all = fill_all(poly, fill_face_ids, fill_face_sizes);
    fill_edges = fill_edges(poly, fill_all);
    poly_wire_faces_impl(poly, pd, bw, s, poly, fid, fill_all); 
    poly_wire_edges_impl(poly, pd, bw, s, poly, fid, fill_edges = fill_edges);
    // dual
    dual_fill_all = fill_all(dual, [], dual_fill_face_sizes);
    dual_fill_edges = fill_edges(dual, dual_fill_all);
    poly_wire_faces_impl(dual, pd, dw0, dscale, poly, fid, dual_fill_all);
    poly_wire_edges_impl(dual, pd, dw0, dscale, poly, fid, fill_edges = dual_fill_edges);
    // wires between primary and dual
    if (cw != 0) {
        cs = (sh - cw) / pd;
        for (i = [0:len(fs) - 1]) {
            c0 = dvs[i] * cs * dfactor;
            for (u = fs[i]) {
                p0 = vs[u] * cs;
                hull() {
                    translate(c0) poly_fill0(poly, cw / pd, fid);
                    translate(p0) poly_fill0(poly, cw / pd, fid);
                }       
            }
        }
    }
}

module poly_wire_dual(
    poly, 
    sh = sh, bw = bw, cw = cw, dw = undef, fid = undef,
    fill_face_ids = [], fill_face_sizes = [], dual_fill_face_sizes = []
) {
    validate(poly);
    pd = diameter(poly, fid);
    translate([0, 0, sh * face_dist(poly, fid) / pd])
        face_rotate(poly, fid) 
            poly_wire_dual0(
                poly, pd, sh, bw, cw, dw, fid, 
                fill_face_ids, fill_face_sizes, dual_fill_face_sizes
            );
}

module validate(poly) {
    vs = poly[0];
    fs = poly[1];
    face_sizes = distinct([
        for (f = fs) len(f)
    ]);
    es = [
        for (f = fs) 
            for (i = [0:len(f) - 1])
                let(j = (i + 1) % len(f))
                    if (f[i] < f[j])
                        [f[i], f[j]]
    ];
    edge_lengths = distinct([
        for (e = es) norm(vs[e[0]] - vs[e[1]])
    ]);
    echo("# Vertices:", N = len(vs));
    echo("# Edges:", N = len(es), edge_lengths = edge_lengths);
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

// octahedron = rectified(tetrahedron);
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

// cube = dual(octahedron);
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

// icosahedron = snub(tetrahedron);
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

// dodecahedron = dual(icosahedron);
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

// --------------------- 13 Arhimedean Solids ---------------------

truncated_tetrahedron = truncated(tetrahedron);
cuboctahedron = rectified(cube);
truncated_cube = truncated(cube);
truncated_octahedron = truncated(octahedron);
rhombicuboctahedron = cantellated(cube);
rhombitruncated_cuboctahedron = beveled(cube);
snub_cube = snub(cube);
icosidodecahedron = rectified(dodecahedron);
truncated_dodecahedron = truncated(dodecahedron);
truncated_icosahedron = truncated(icosahedron);
rhombicosidodecahedron = cantellated(dodecahedron);
rhombitruncated_icosidodecahedron = beveled(dodecahedron);
snub_dodecahedron = snub(dodecahedron);

// --------------------- 13 Catalan Solids (Arhimedean Duals) ---------------------

triakis_tetrahedron = dual(truncated_tetrahedron);
rhombic_dodecahedron = dual(cuboctahedron);
triakis_octahedron = dual(truncated_cube);
tetrakis_hexahedron = dual(truncated_octahedron);
deltoidal_icositetrahedron = dual(rhombicuboctahedron);
disdyakis_dodecahedron = dual(rhombitruncated_cuboctahedron);
pentagonal_icositetrahedron = dual(snub_cube);
rhombic_triacontahedron = dual(icosidodecahedron);
triakis_icosahedron = dual(truncated_dodecahedron);
pentakis_dodecahedron = dual(truncated_icosahedron);
deltoidal_hexecontahedron = dual(rhombicosidodecahedron);
disdyakis_triacontahedron = dual(rhombitruncated_icosidodecahedron);
pentagonal_hexecontahedron = dual(snub_dodecahedron);

// --------------------- Stellated polyhedra ---------------------

star_equilateral_24deltahedron= augmented(cube, 1/sqrt(2));
star_equilateral_60deltahedron = augmented(dodecahedron, sqrt(0.5 - 0.1 * sqrt(5)));
small_stellated_dodecahedron = augmented(dodecahedron, sqrt(1 + 2/5 * sqrt(5)));
small_triambic_icosahedron = augmented(icosahedron, sqrt(15)/15);
star_equilateral_60deltahedron2 = augmented(icosahedron, sqrt(6)/3);
great_stellated_dodecahedron = augmented(icosahedron, sqrt(3)*(3 + sqrt(5))/6);
