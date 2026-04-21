pub fn shapefunc<Real>(pco: &[Real; 3]) -> [Real; 6]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let r = pco[0];
    let s = pco[1];
    let t = pco[2];
    let tm = one - t;
    let rs = one - r - s;
    [rs * tm, r * tm, s * tm, rs * t, r * t, s * t]
}

/// Returns parametric coordinates (r,s,t) of the origin inside the prism, or None if outside.
/// Valid domain: r>=0, s>=0, r+s<=1, t in [0,1].
pub fn parametric_coord_for_origin<Real>(
    p0: &[Real; 3], p1: &[Real; 3], p2: &[Real; 3],
    p3: &[Real; 3], p4: &[Real; 3], p5: &[Real; 3],
) -> Option<[Real; 3]>
where
    Real: num_traits::Float,
{
    use slice_of_array::SliceFlatExt;
    let zero = Real::zero();
    let one = Real::one();
    let third = one / (one + one + one);
    let half = one / (one + one);
    let mut pco = [third, third, half];
    for _ in 0..20 {
        let sf = shapefunc(&pco);
        let pos = [
            sf[0]*p0[0] + sf[1]*p1[0] + sf[2]*p2[0] + sf[3]*p3[0] + sf[4]*p4[0] + sf[5]*p5[0],
            sf[0]*p0[1] + sf[1]*p1[1] + sf[2]*p2[1] + sf[3]*p3[1] + sf[4]*p4[1] + sf[5]*p5[1],
            sf[0]*p0[2] + sf[1]*p1[2] + sf[2]*p2[2] + sf[3]*p3[2] + sf[4]*p4[2] + sf[5]*p5[2],
        ];
        let jac = dxdr(p0, p1, p2, p3, p4, p5, &pco);
        let jac_col: &[Real; 9] = jac.flat().try_into().unwrap();
        let jac_col = crate::mat3_col_major::transpose(jac_col);
        let j_inv = crate::mat3_col_major::try_inverse(&jac_col)?;
        let dpco = crate::mat3_col_major::mult_vec(&j_inv, &pos);
        for i in 0..3 { pco[i] = pco[i] - dpco[i]; }
    }
    let (r, s, t) = (pco[0], pco[1], pco[2]);
    if r >= zero && s >= zero && r + s <= one && t >= zero && t <= one {
        Some(pco)
    } else {
        None
    }
}

/// Returns the nearest point to the origin on the prism surface, and the 6 shape function weights.
/// Vertices: p0,p1,p2 are the bottom triangle, p3,p4,p5 are the top triangle.
pub fn nearest_to_origin<Real>(
    p0: &[Real; 3], p1: &[Real; 3], p2: &[Real; 3],
    p3: &[Real; 3], p4: &[Real; 3], p5: &[Real; 3],
) -> ([Real; 3], [Real; 6])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let zero = Real::zero();
    let one = Real::one();

    if let Some(pco) = parametric_coord_for_origin(p0, p1, p2, p3, p4, p5) {
        return ([zero; 3], shapefunc(&pco));
    }

    let mut p_min = *p0;
    let mut w_min = [one, zero, zero, zero, zero, zero];
    let mut d_min = crate::vec3::norm(p0);

    let mut update = |p: [Real; 3], w: [Real; 6]| {
        let d = crate::vec3::norm(&p);
        if d < d_min {
            d_min = d;
            p_min = p;
            w_min = w;
        }
    };

    // triangular faces
    {
        let (p, a0, a1, a2) = crate::tri3::nearest_to_origin3(p0, p1, p2);
        update(p, [a0, a1, a2, zero, zero, zero]);
    }
    {
        let (p, a3, a4, a5) = crate::tri3::nearest_to_origin3(p3, p4, p5);
        update(p, [zero, zero, zero, a3, a4, a5]);
    }

    // quad faces: weights follow quad::position formula w0=(1-s0)(1-s1), w1=s0(1-s1), w2=s0*s1, w3=(1-s0)*s1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p0, p1, p4, p3);
        let (w0, w1, w2, w3) = ((one-s0)*(one-s1), s0*(one-s1), s0*s1, (one-s0)*s1);
        update(p, [w0, w1, zero, w3, w2, zero]); // q0→p0, q1→p1, q2→p4, q3→p3
    }
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p1, p2, p5, p4);
        let (w0, w1, w2, w3) = ((one-s0)*(one-s1), s0*(one-s1), s0*s1, (one-s0)*s1);
        update(p, [zero, w0, w1, zero, w3, w2]); // q0→p1, q1→p2, q2→p5, q3→p4
    }
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p2, p0, p3, p5);
        let (w0, w1, w2, w3) = ((one-s0)*(one-s1), s0*(one-s1), s0*s1, (one-s0)*s1);
        update(p, [w1, zero, w0, w2, zero, w3]); // q0→p2, q1→p0, q2→p3, q3→p5
    }

    // edges
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p0, p1); update(p, [s0, s1, zero, zero, zero, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p1, p2); update(p, [zero, s0, s1, zero, zero, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p2, p0); update(p, [s1, zero, s0, zero, zero, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p3, p4); update(p, [zero, zero, zero, s0, s1, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p4, p5); update(p, [zero, zero, zero, zero, s0, s1]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p5, p3); update(p, [zero, zero, zero, s1, zero, s0]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p0, p3); update(p, [s0, zero, zero, s1, zero, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p1, p4); update(p, [zero, s0, zero, zero, s1, zero]); }
    { let (p, s0, s1) = crate::edge3::nearest_to_origin3(p2, p5); update(p, [zero, zero, s0, zero, zero, s1]); }

    (p_min, w_min)
}

#[test]
fn test_shapefunc() {
    let _a = shapefunc(&[0., 0., 0.]);
    let _a = shapefunc(&[1., 0., 0.]);
    let _a = shapefunc(&[0., 1., 0.]);
    let _a = shapefunc(&[0., 0., 1.]);
}

pub fn dndr<Real>(pco: &[Real; 3]) -> [[Real; 3]; 6]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let r = pco[0];
    let s = pco[1];
    let t = pco[2];
    [
        [-one + t, -one + t, -one + r + s],
        [one - t, zero, -r],
        [zero, one - t, -s],
        [-t, -t, one - r - s],
        [t, zero, r],
        [zero, t, s],
    ]
}

pub fn dxdr<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    pco: &[Real; 3],
) -> [[Real; 3]; 3]
where
    Real: num_traits::Float,
{
    let dndr = dndr(pco);
    let mut dxdr = [[Real::zero(); 3]; 3];
    for idim in 0..3 {
        for ir in 0..3 {
            dxdr[idim][ir] = dxdr[idim][ir] + p0[idim] * dndr[0][ir];
            dxdr[idim][ir] = dxdr[idim][ir] + p1[idim] * dndr[1][ir];
            dxdr[idim][ir] = dxdr[idim][ir] + p2[idim] * dndr[2][ir];
            dxdr[idim][ir] = dxdr[idim][ir] + p3[idim] * dndr[3][ir];
            dxdr[idim][ir] = dxdr[idim][ir] + p4[idim] * dndr[4][ir];
            dxdr[idim][ir] = dxdr[idim][ir] + p5[idim] * dndr[5][ir];
        }
    }
    dxdr
}

#[test]
fn test_dxdr() {
    let p0 = [0., 0., 0.];
    let p1 = [1., 0., 0.];
    let p2 = [0., 1., 0.];
    let p3 = [0., 0., 1.];
    let p4 = [1., 0., 1.];
    let p5 = [0., 1., 1.];
    let dxdr: [[f64; 3]; 3] = dxdr(&p0, &p1, &p2, &p3, &p4, &p5, &[1. / 3., 1. / 3., 0.5]);
    use slice_of_array::SliceFlatExt;
    let vol = crate::mat3_col_major::determinant(dxdr.flat().try_into().unwrap());
    assert!((vol - 1.0).abs() < 1.0e-10);
}

pub fn volume<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    i_gauss_degree: usize,
) -> Real
where
    Real: num_traits::Float + 'static + std::fmt::Debug,
    crate::quadrature_line::Quad<Real>: crate::quadrature_line::QuadratureLine<Real>,
    crate::quadrature_tri::Quad<Real>: crate::quadrature_tri::QuadratureTri<Real>,
{
    let one = Real::one();
    let half = one / (one + one);
    let one4th = half * half;
    use crate::quadrature_line::QuadratureLine;
    use crate::quadrature_tri::QuadratureTri;
    let quad_l: &[[Real; 2]] = crate::quadrature_line::Quad::<Real>::hoge(i_gauss_degree);
    let quad_t: &[[Real; 3]] = crate::quadrature_tri::Quad::<Real>::hoge(i_gauss_degree);
    let mut volume = Real::zero();
    let num_quad_l = quad_l.len();
    let num_quad_t = quad_t.len();
    for ir_l in 0..num_quad_l {
        for ir_t in 0..num_quad_t {
            let pco = [
                quad_t[ir_t][0],
                quad_t[ir_t][1],
                (quad_l[ir_l][0] + one) * half,
            ];
            let w = quad_l[ir_l][1] * quad_t[ir_t][2];
            let dxdr = dxdr(p0, p1, p2, p3, p4, p5, &pco);
            use slice_of_array::SliceFlatExt;
            let detjac =
                crate::mat3_col_major::determinant::<Real>(dxdr.flat().try_into().unwrap());
            volume = volume + detjac * w;
        }
    }
    volume * one4th
}

#[test]
fn test_nearest_to_origin() {
    // prism with bottom triangle in z=2 plane, top in z=4 plane
    // the origin projects onto the interior of the bottom face: (0,0,2) has bary coords (0.25,0.25,0.5)
    let p0 = [-1.0f64, -1.0, 2.0];
    let p1 = [ 1.0, -1.0, 2.0];
    let p2 = [ 0.0,  1.0, 2.0];
    let p3 = [-1.0, -1.0, 4.0];
    let p4 = [ 1.0, -1.0, 4.0];
    let p5 = [ 0.0,  1.0, 4.0];

    {
        let v = volume(&p0, &p1, &p2, &p3, &p4, &p5, 2);
        dbg!(v);
    }

    let (p, w) = nearest_to_origin(&p0, &p1, &p2, &p3, &p4, &p5);

    // nearest point should be on the bottom face at (0, 0, 2)
    assert!((p[0]).abs() < 1.0e-10 && (p[1]).abs() < 1.0e-10 && (p[2] - 2.0).abs() < 1.0e-10,
        "expected (0,0,2), got {:?}", p);

    // weights must sum to 1 and be non-negative
    let wsum: f64 = w.iter().sum();
    assert!((wsum - 1.0).abs() < 1.0e-10, "weights sum to {}", wsum);
    assert!(w.iter().all(|&wi| wi >= -1.0e-10), "negative weight: {:?}", w);

    // reconstruction: sum(w[i] * p_i) == p
    let verts = [p0, p1, p2, p3, p4, p5];
    let recon: [f64; 3] = std::array::from_fn(|i| w.iter().zip(verts.iter()).map(|(&wi, vi)| wi * vi[i]).sum());
    assert!(crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10, "reconstruction failed");

    // nearest point must be no farther than any vertex
    let d = crate::vec3::norm(&p);
    for v in &verts {
        assert!(d <= crate::vec3::norm(v) + 1.0e-10, "vertex {:?} is closer than result", v);
    }
}

