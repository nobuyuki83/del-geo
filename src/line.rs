//! methods for 2D or 3D line

use num_traits::AsPrimitive;

pub fn nearest_to_point<T, const X: usize>(
    p: &nalgebra::base::SVector<T, X>, // point
    s: &nalgebra::base::SVector<T, X>, // source
    d: &nalgebra::base::SVector<T, X>,
) -> (nalgebra::base::SVector<T, X>, T)
where
    T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>,
{
    if d.dot(d) < 1.0e-20_f64.as_() {
        return (*s, 0_f64.as_());
    }
    let a = d.dot(d);
    let b = d.dot(&(s - p));
    let t: T = -b / a;
    (s + d.scale(t), t)
}
