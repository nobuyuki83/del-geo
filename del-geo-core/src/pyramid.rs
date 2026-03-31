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
