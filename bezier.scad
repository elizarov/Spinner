//$fs = $preview ? 1.0 : 0.1;

//bezier_surface([
//    [ [0, -10,  0], [10, 0,  0], [0, 10,  0], [-10, 0,  0] ],
//    [ [0, -10, 10], [ 5, 0, 10], [8, 20, 10], [-20, 0, 10] ],
//    [ [0, -10, 20], [10, 0, 20], [0, 10, 20], [-10, 0, 20] ],
//]);


// quadratic Bézier curves
module bezier_surface(p) {
    assert(len(p) % 2 == 1, "Expected odd number of points");
    n = (len(p) - 1) / 2;
    k = len(p[0]);
    steps = [
        for (i = [0:n-1]) 
            max([
                for (j = [0:k-1]) 
                let(p0 = p[2 * i][j])
                let(p1 = p[2 * i + 1][j])
                let(p2 = p[2 * i + 2][j])
                let(d = norm(p2 - p1) + norm(p1 - p0))
                ceil(d / $fs)
            ])
    ];
    v = [
        for (i = [0:n-1]) 
            let(st = steps[i])
            for (ti = (i == 0 ? [0:st] : [1:st]))
                let(t = ti / st)
                let(t1 = 1 - t)
                for (j = [0:k-1]) 
                    let(p0 = p[2 * i][j])
                    let(p1 = p[2 * i + 1][j])
                    let(p2 = p[2 * i + 2][j])
                        t1 * t1 * p0 + 2 * t1 * t * p1 + t * t * p2
    ];
    m = len(v);            
    f0 = [for(j = [0:k-1]) j];
    f1 = [for(j = [m-1:-1:m - k]) j];
    f2 = [
        for (i = [1:m/k-1]) 
            for(j = [0:k-1])
            [ 
                i * k + j, 
                i * k + (j + 1) % k, 
                (i - 1) * k + (j + 1) % k, 
                (i - 1) * k + j
            ]           
    ];
    polyhedron(v, concat([f0, f1], f2));  
}