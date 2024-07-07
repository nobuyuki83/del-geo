//! functions for N-dimensional edge (line segment)

use num_traits::AsPrimitive;

pub fn nearest_to_origin<T, const X: usize>(
    p0: &nalgebra::base::SVector<T, X>, // start
    p1: &nalgebra::base::SVector<T, X>,
) -> nalgebra::base::SVector<T, X>
where
    T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>,
{
    let d = p1 - p0;
    let a: T = d.dot(&d);
    if a < 1.0e-20_f64.as_() {
        return (p0 + p1) * 0.5_f64.as_();
    }
    let b = d.dot(p0);
    let mut r0: T = -b / a;
    if r0 < T::zero() {
        r0 = T::zero();
    }
    if r0 > T::one() {
        r0 = T::one();
    }
    p0.scale(T::one() - r0) + p1.scale(r0)
}

pub fn distance_to_point<T, const X: usize>(
    po_c: &nalgebra::base::SVector<T, X>,
    po_s: &nalgebra::base::SVector<T, X>,
    po_e: &nalgebra::base::SVector<T, X>,
) -> T
where
    T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>,
{
    crate::edge::nearest_to_origin(&(po_s - po_c), &(po_e - po_c)).norm()
}
