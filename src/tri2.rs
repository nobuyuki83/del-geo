//! methods for 2D triangle

use num_traits::AsPrimitive;

pub fn area_<T>(
    p0: &[T; 2],
    p1: &[T; 2],
    p2: &[T; 2]) -> T
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + 'static + Copy,
          f64: AsPrimitive<T>
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    0.5_f64.as_() * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
}

/// shape function's derivative in x- and y-direction and the constant term
///
/// # Example
/// ```
/// let (dldx, c) = del_geo::tri2::dldx_(&[0., 1.], &[1., 2.], &[3., 4.]);
/// let a = dldx[0][1]; // x-differentiation of L^1
/// ```
pub fn dldx_<T>(
    p0: &[T; 2],
    p1: &[T; 2],
    p2: &[T; 2]) -> ([[T; 3]; 2], [T; 3])
    where T: num_traits::Float + Copy + 'static,
          f64: AsPrimitive<T>
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let a0 = area_(p0, p1, p2);
    let tmp1: T = 0.5.as_() / a0;
    (
        [
            [
                tmp1 * (p1[1] - p2[1]),
                tmp1 * (p2[1] - p0[1]),
                tmp1 * (p0[1] - p1[1])],
            [
                tmp1 * (p2[0] - p1[0]),
                tmp1 * (p0[0] - p2[0]),
                tmp1 * (p1[0] - p0[0])]
        ],
        [
            tmp1 * (p1[0] * p2[1] - p2[0] * p1[1]),
            tmp1 * (p2[0] * p0[1] - p0[0] * p2[1]),
            tmp1 * (p0[0] * p1[1] - p1[0] * p0[1])]
    )
}


