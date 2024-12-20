//! methods for 2D triangle

pub fn area<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    half * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
}

/// # Return
/// `(dldp0: [T;2], dldp1: [T;2], dldp2: [T;2])`
pub fn dldw_area<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2], dldarea: T) -> ([T; 2], [T; 2], [T; 2])
where
    T: num_traits::Float,
{
    let two = T::one() + T::one();
    let half = T::one() / two;
    let dareadp0x2 = [p1[1] - p2[1], p2[0] - p1[0]];
    let dareadp1x2 = [p2[1] - p0[1], p0[0] - p2[0]];
    let dareadp2x2 = [p0[1] - p1[1], p1[0] - p0[0]];
    (
        crate::vec2::scale(&dareadp0x2, dldarea * half),
        crate::vec2::scale(&dareadp1x2, dldarea * half),
        crate::vec2::scale(&dareadp2x2, dldarea * half),
    )
}

#[test]
fn test_dldw_area() {
    let p0s = [[0.1, -0.2], [1.3, 0.2], [0.6, 0.45]];
    let dldarea = 1.3f32;
    let area0 = area(&p0s[0], &p0s[1], &p0s[2]);
    let l0 = dldarea * area0;
    let dl = dldw_area(&p0s[0], &p0s[1], &p0s[2], dldarea);
    let dl = [dl.0, dl.1, dl.2];
    let eps = 1.0e-4;
    for inode in 0..3 {
        for idim in 0..2 {
            let p1s = {
                let mut p1s = p0s;
                p1s[inode][idim] += eps;
                p1s
            };
            let area1 = area(&p1s[0], &p1s[1], &p1s[2]);
            let l1 = dldarea * area1;
            let val_num = (l1 - l0) / eps;
            let val_ana = dl[inode][idim];
            let diff_abs = (val_num - val_ana).abs();
            assert!(diff_abs < 5.0e-4, "{}", diff_abs);
        }
    }
}

pub fn is_inside<Real>(
    p0: &[Real; 2],
    p1: &[Real; 2],
    p2: &[Real; 2],
    q: &[Real; 2],
    sign: Real,
) -> Option<(Real, Real)>
where
    Real: num_traits::Float,
{
    let a0 = area(q, p1, p2) * sign;
    if a0 < Real::zero() {
        return None;
    }
    let a1 = area(q, p2, p0) * sign;
    if a1 < Real::zero() {
        return None;
    }
    let a2 = area(q, p0, p1) * sign;
    if a2 < Real::zero() {
        return None;
    }
    let sum_area_inv = Real::one() / (a0 + a1 + a2);
    Some((a0 * sum_area_inv, a1 * sum_area_inv))
}

/// shape function's derivative in x- and y-direction and the constant term
///
/// # Example
/// ```
/// let (dldx, c) = del_geo_core::tri2::dldx(&[0., 1.], &[1., 2.], &[3., 4.]);
/// let a = dldx[0][2]; // x-differentiation of L^2
/// ```
pub fn dldx<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> ([[T; 3]; 2], [T; 3])
where
    T: num_traits::Float,
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let half = T::one() / (T::one() + T::one());
    let a0 = area(p0, p1, p2);
    let tmp1: T = half / a0;
    (
        [
            [
                tmp1 * (p1[1] - p2[1]),
                tmp1 * (p2[1] - p0[1]),
                tmp1 * (p0[1] - p1[1]),
            ],
            [
                tmp1 * (p2[0] - p1[0]),
                tmp1 * (p0[0] - p2[0]),
                tmp1 * (p1[0] - p0[0]),
            ],
        ],
        [
            tmp1 * (p1[0] * p2[1] - p2[0] * p1[1]),
            tmp1 * (p2[0] * p0[1] - p0[0] * p2[1]),
            tmp1 * (p0[0] * p1[1] - p1[0] * p0[1]),
        ],
    )
}

pub fn barycentric_coords<Real>(
    p0: &[Real; 2],
    p1: &[Real; 2],
    p2: &[Real; 2],
    q: &[Real; 2],
) -> Option<(Real, Real, Real)>
where
    Real: num_traits::Float,
{
    let a0 = area(q, p1, p2);
    let a1 = area(q, p2, p0);
    let a2 = area(q, p0, p1);
    if (a0 + a1 + a2).is_zero() {
        return None;
    }
    let sum_area_inv = Real::one() / (a0 + a1 + a2);
    Some((a0 * sum_area_inv, a1 * sum_area_inv, a2 * sum_area_inv))
}

// -------------------------------------------

pub struct Tri2<'a, Real> {
    pub p0: &'a [Real; 2],
    pub p1: &'a [Real; 2],
    pub p2: &'a [Real; 2],
}

#[allow(clippy::needless_lifetimes)]
impl<'a, Real> Tri2<'a, Real>
where
    Real: num_traits::Float,
{
    pub fn is_inside(&self, q: &[Real; 2], sign: Real) -> Option<(Real, Real)> {
        is_inside(self.p0, self.p1, self.p2, q, sign)
    }

    pub fn area(&self) -> Real {
        area(self.p0, self.p1, self.p2)
    }
}
