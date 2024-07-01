//! methods for 2D triangle

pub fn area<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let half = T::one() / (T::one() + T::one());
    half * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
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
/// let (dldx, c) = del_geo_core::tri2::dldx_(&[0., 1.], &[1., 2.], &[3., 4.]);
/// let a = dldx[0][2]; // x-differentiation of L^2
/// ```
pub fn dldx_<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> ([[T; 3]; 2], [T; 3])
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
