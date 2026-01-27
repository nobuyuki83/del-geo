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
    [rs * tm, s * tm, r * tm, rs * t, s * t, r * t]
}

#[test]
fn test_shapefunc() {
    let a = shapefunc(&[0., 0., 0.]);
    let a = shapefunc(&[1., 0., 0.]);
    let a = shapefunc(&[0., 1., 0.]);
    let a = shapefunc(&[0., 0., 1.]);
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
        [zero, one - t, -s],
        [one - t, zero, -r],
        [-t, -t, one - r - s],
        [zero, t, s],
        [t, zero, r],
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
    let p1 = [0., 1., 0.];
    let p2 = [1., 0., 0.];
    let p3 = [0., 0., 1.];
    let p4 = [0., 1., 1.];
    let p5 = [1., 0., 1.];
    let dxdr = dxdr(&p0, &p1, &p2, &p3, &p4, &p5, &[1. / 3., 1. / 3., 0.5]);
    dbg!(dxdr);
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
