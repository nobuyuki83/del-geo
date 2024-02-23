//! methods for 2D triangle

use num_traits::AsPrimitive;

pub fn area_<T>(
    p0: &[T; 2],
    p1: &[T; 2],
    p2: &[T; 2]) -> T
    where T: num_traits::Float
{
    assert!(p0.len() == 2 && p1.len() == 2 && p2.len() == 2);
    let half = T::one() / (T::one() + T::one());
    half * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
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

//////////


pub fn area<T>(
    v1: &nalgebra::Vector2<T>,
    v2: &nalgebra::Vector2<T>,
    v3: &nalgebra::Vector2<T>) -> T
    where T: num_traits::Float + Copy,
{
    ((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])) / (T::one()+T::one())
}

pub fn circumcenter<T>(
    p0: &nalgebra::Vector2<T>,
    p1: &nalgebra::Vector2<T>,
    p2: &nalgebra::Vector2<T>) -> nalgebra::Vector2::<T>
where T: num_traits::Float + Copy + std::fmt::Debug,
{
    let area = crate::tri2::area(p0, p1, p2);
    if area.is_zero() { panic!(); }
    let four = T::one() + T::one() + T::one() + T::one();
    let sixteen = four * four;
    let tmp_val = T::one() / (area * area * sixteen);

    let dtmp0 = crate::edge2::length_squared(p1,p2);
    let dtmp1 = crate::edge2::length_squared(p0,p2);
    let dtmp2 = crate::edge2::length_squared(p0,p1);

    let etmp0: T = tmp_val * dtmp0 * (dtmp1 + dtmp2 - dtmp0);
    let etmp1: T = tmp_val * dtmp1 * (dtmp0 + dtmp2 - dtmp1);
    let etmp2: T = tmp_val * dtmp2 * (dtmp0 + dtmp1 - dtmp2);

    nalgebra::Vector2::<T>::new(
        etmp0 * p0[0] + etmp1 * p1[0] + etmp2 * p2[0],
        etmp0 * p0[1] + etmp1 * p1[1] + etmp2 * p2[1])
}