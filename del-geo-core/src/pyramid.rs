//! vtx index is the same as the VTK

pub fn shapefunc<Real>(pco: [Real; 3]) -> [Real; 5]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let r = pco[0];
    let s = pco[1];
    let t = pco[2];

    let rm = one - r;
    let sm = one - s;
    let tm = one - t;

    [
        rm * sm * tm, // sf[0]
        r * sm * tm,  // sf[1]
        r * s * tm,   // sf[2]
        rm * s * tm,  // sf[3]
        t,            // sf[4]
    ]
}

/// Returns parametric coordinates of the origin inside the pyramid, or None if outside.
fn parametric_coord_for_origin<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
) -> Option<[Real; 3]>
where
    Real: num_traits::Float,
{
    use slice_of_array::SliceFlatExt;
    let zero = Real::zero();
    let one = Real::one();
    let half = one / (one + one);
    let mut pco = [half; 3];
    for _ in 0..20 {
        let sf = shapefunc(pco);
        let pos = [
            sf[0] * p0[0] + sf[1] * p1[0] + sf[2] * p2[0] + sf[3] * p3[0] + sf[4] * p4[0],
            sf[0] * p0[1] + sf[1] * p1[1] + sf[2] * p2[1] + sf[3] * p3[1] + sf[4] * p4[1],
            sf[0] * p0[2] + sf[1] * p1[2] + sf[2] * p2[2] + sf[3] * p3[2] + sf[4] * p4[2],
        ];
        let jac = dxdr(p0, p1, p2, p3, p4, &pco);
        let jac_col: &[Real; 9] = jac.flat().try_into().unwrap();
        let jac_col = crate::mat3_col_major::transpose(jac_col);
        let j_inv = crate::mat3_col_major::try_inverse(&jac_col)?;
        let dpco = crate::mat3_col_major::mult_vec(&j_inv, &pos);
        for i in 0..3 {
            pco[i] = pco[i] - dpco[i];
        }
    }
    if pco[0] >= zero
        && pco[0] <= one
        && pco[1] >= zero
        && pco[1] <= one
        && pco[2] >= zero
        && pco[2] <= one
    {
        Some(pco)
    } else {
        None
    }
}

/// Returns the nearest point to the origin on the pyramid surface, and the 5 shape function weights.
/// Vertices: p0,p1,p2,p3 are the base quad (VTK order), p4 is the apex.
pub fn nearest_to_origin<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
) -> ([Real; 3], [Real; 5])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let zero = Real::zero();
    let one = Real::one();

    // if origin is inside the pyramid, return origin with interior shape function weights
    if let Some(pco) = parametric_coord_for_origin(p0, p1, p2, p3, p4) {
        let sf = shapefunc(pco);
        return ([zero; 3], sf);
    }

    let mut p_min = *p0;
    let mut w_min = [one, zero, zero, zero, zero];
    let mut d_min = crate::vec3::norm(p0);

    let mut update = |p: [Real; 3], w: [Real; 5]| {
        let d = crate::vec3::norm(&p);
        if d < d_min {
            d_min = d;
            p_min = p;
            w_min = w;
        }
    };

    // base quad face: w0=(1-s0)(1-s1)→p0, w1=s0(1-s1)→p1, w2=s0*s1→p2, w3=(1-s0)*s1→p3
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p0, p1, p2, p3);
        let (w0, w1, w2, w3) = (
            (one - s0) * (one - s1),
            s0 * (one - s1),
            s0 * s1,
            (one - s0) * s1,
        );
        update(p, [w0, w1, w2, w3, zero]);
    }

    // triangular lateral faces
    {
        let (p, a0, a1, a4) = crate::tri3::nearest_to_origin3(p0, p1, p4);
        update(p, [a0, a1, zero, zero, a4]);
    }
    {
        let (p, a1, a2, a4) = crate::tri3::nearest_to_origin3(p1, p2, p4);
        update(p, [zero, a1, a2, zero, a4]);
    }
    {
        let (p, a2, a3, a4) = crate::tri3::nearest_to_origin3(p2, p3, p4);
        update(p, [zero, zero, a2, a3, a4]);
    }
    {
        let (p, a3, a0, a4) = crate::tri3::nearest_to_origin3(p3, p0, p4);
        update(p, [a0, zero, zero, a3, a4]);
    }

    // base edges
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p0, p1);
        update(p, [s0, s1, zero, zero, zero]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p1, p2);
        update(p, [zero, s0, s1, zero, zero]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p2, p3);
        update(p, [zero, zero, s0, s1, zero]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p3, p0);
        update(p, [s1, zero, zero, s0, zero]);
    }

    // lateral edges
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p0, p4);
        update(p, [s0, zero, zero, zero, s1]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p1, p4);
        update(p, [zero, s0, zero, zero, s1]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p2, p4);
        update(p, [zero, zero, s0, zero, s1]);
    }
    {
        let (p, s0, s1) = crate::edge3::nearest_to_origin3(p3, p4);
        update(p, [zero, zero, zero, s0, s1]);
    }

    (p_min, w_min)
}

/// pco: parametric coordinates
pub fn dndr<Real>(pco: &[Real; 3]) -> [[Real; 3]; 5]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let r = pco[0];
    let s = pco[1];
    let t = pco[2];

    let rm = one - r;
    let sm = one - s;
    let tm = one - t;

    let mut dndr = [[zero; 3]; 5];

    // r-derivatives
    dndr[0][0] = -sm * tm;
    dndr[1][0] = sm * tm;
    dndr[2][0] = s * tm;
    dndr[3][0] = -s * tm;
    dndr[4][0] = zero;

    // s-derivatives
    dndr[0][1] = -rm * tm;
    dndr[1][1] = -r * tm;
    dndr[2][1] = r * tm;
    dndr[3][1] = rm * tm;
    dndr[4][1] = zero;

    // t-derivatives
    dndr[0][2] = -rm * sm;
    dndr[1][2] = -r * sm;
    dndr[2][2] = -r * s;
    dndr[3][2] = -rm * s;
    dndr[4][2] = one;

    dndr
}

pub fn dxdr<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
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
        }
    }
    dxdr
}

pub fn volume<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    i_gauss_degree: usize,
) -> Real
where
    Real: num_traits::Float + 'static + std::fmt::Debug,
    crate::quadrature_line::Quad<Real>: crate::quadrature_line::QuadratureLine<Real>,
{
    let one = Real::one();
    let half = one / (one + one);
    let one8th = half * half * half;
    use crate::quadrature_line::QuadratureLine;
    let quadrature: &[[Real; 2]] = crate::quadrature_line::Quad::<Real>::hoge(i_gauss_degree);
    let mut volume = Real::zero();
    let num_quadr = quadrature.len();
    for (ir1, ir2, ir3) in itertools::iproduct!(0..num_quadr, 0..num_quadr, 0..num_quadr) {
        let pco = [
            (quadrature[ir1][0] + one) * half,
            (quadrature[ir2][0] + one) * half,
            (quadrature[ir3][0] + one) * half,
        ];
        let w = quadrature[ir1][1] * quadrature[ir2][1] * quadrature[ir3][1];
        let dxdr = dxdr(p0, p1, p2, p3, p4, &pco);
        use slice_of_array::SliceFlatExt;
        let detjac = crate::mat3_col_major::determinant::<Real>(dxdr.flat().try_into().unwrap());
        volume = volume + detjac * w;
    }
    volume * one8th
}

#[test]
fn test_nearest_to_origin() {
    // origin inside: base quad at z=-1, apex at z=2
    {
        let p0 = [-1.0f64, -1.0, -1.0];
        let p1 = [1.0, -1.0, -1.0];
        let p2 = [1.0, 1.0, -1.0];
        let p3 = [-1.0, 1.0, -1.0];
        let p4 = [0.0, 0.0, 2.0];
        let (p, w) = nearest_to_origin(&p0, &p1, &p2, &p3, &p4);
        assert!(
            crate::vec3::norm(&p) < 1.0e-10,
            "inside: nearest point should be origin, got {:?}",
            p
        );
        assert!(
            w.iter().all(|&wi| wi > 0.0),
            "inside: all weights should be positive: {:?}",
            w
        );
        let wsum: f64 = w.iter().sum();
        assert!((wsum - 1.0).abs() < 1.0e-10, "weights sum to {}", wsum);
        let vs = [p0, p1, p2, p3, p4];
        let recon: [f64; 3] =
            std::array::from_fn(|i| w.iter().zip(vs.iter()).map(|(&wi, vi)| wi * vi[i]).sum());
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction failed"
        );
    }

    // origin outside: pyramid shifted to x>=3, nearest point on the base face
    {
        let p0 = [3.0f64, -1.0, -1.0];
        let p1 = [5.0, -1.0, -1.0];
        let p2 = [5.0, 1.0, -1.0];
        let p3 = [3.0, 1.0, -1.0];
        let p4 = [4.0, 0.0, 2.0];
        let (p, w) = nearest_to_origin(&p0, &p1, &p2, &p3, &p4);
        let wsum: f64 = w.iter().sum();
        assert!((wsum - 1.0).abs() < 1.0e-10, "weights sum to {}", wsum);
        assert!(
            w.iter().all(|&wi| wi >= -1.0e-10),
            "negative weight: {:?}",
            w
        );
        let vs = [p0, p1, p2, p3, p4];
        let recon: [f64; 3] =
            std::array::from_fn(|i| w.iter().zip(vs.iter()).map(|(&wi, vi)| wi * vi[i]).sum());
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction failed"
        );
        let d = crate::vec3::norm(&p);
        for v in &vs {
            assert!(
                d <= crate::vec3::norm(v) + 1.0e-10,
                "vertex {:?} is closer than result",
                v
            );
        }
    }
}
