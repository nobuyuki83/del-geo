//! methods for 3D edge (line segment)

use num_traits::AsPrimitive;
/// trait for 3D edge (line segment)
pub trait Edge3Trait<T> {
    fn length(&self, other: &Self) -> T;
    fn squared_length(&self, other: &Self) -> T;
    fn nearest_to_point3(&self, p1: &Self, point_pos: &Self) -> Self;
}
impl<Real> Edge3Trait<Real> for [Real; 3]
where
    Real: num_traits::Float + 'static + std::ops::MulAssign,
    f64: AsPrimitive<Real>,
{
    fn length(&self, other: &Self) -> Real {
        length(self, other)
    }
    fn squared_length(&self, other: &Self) -> Real {
        squared_length(self, other)
    }
    fn nearest_to_point3(&self, p1: &Self, point_pos: &Self) -> Self {
        nearest_to_point3(self, p1, point_pos)
    }
}
pub fn length<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let x = p0[0] - p1[0];
    let y = p0[1] - p1[1];
    let z = p0[2] - p1[2];
    (x * x + y * y + z * z).sqrt()
}

pub fn squared_length<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let x = p0[0] - p1[0];
    let y = p0[1] - p1[1];
    let z = p0[2] - p1[2];
    x * x + y * y + z * z
}

pub fn nearest_to_point3<T>(p0: &[T; 3], p1: &[T; 3], point_pos: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float + 'static + Copy + PartialOrd + std::ops::MulAssign,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let d = std::array::from_fn(|i| p1[i] - p0[i]);
    let t = {
        if d.dot(&d) > T::epsilon() {
            let ps = std::array::from_fn(|i| p0[i] - point_pos[i]);
            let a = d.dot(&d);
            let b = d.dot(&ps);
            (-b / a).clamp(0f64.as_(), 1f64.as_())
        } else {
            0.5_f64.as_()
        }
    };
    std::array::from_fn(|i| p0[i] + t * d[i])
}
