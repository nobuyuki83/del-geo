//! methods for 3D edge (line segment)

use num_traits::AsPrimitive;

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
    T: num_traits::Float + 'static + Copy + PartialOrd,
    f64: AsPrimitive<T>,
{
    use crate::vec3;
    let d = std::array::from_fn(|i| p1[i] - p0[i]);
    let t = {
        if vec3::dot(&d, &d) > T::epsilon() {
            let ps = std::array::from_fn(|i| p0[i] - point_pos[i]);
            let a = vec3::dot(&d, &d);
            let b = vec3::dot(&d, &ps);
            (-b / a).clamp(0f64.as_(), 1f64.as_())
        } else {
            0.5_f64.as_()
        }
    };
    std::array::from_fn(|i| p0[i] + t * d[i])
}
