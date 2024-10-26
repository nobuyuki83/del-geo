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

pub fn nearest_to_point3<T>(p0: &[T; 3], p1: &[T; 3], point_pos: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float + 'static + Copy + PartialOrd,
    f64: AsPrimitive<T>,
{
    use crate::vec3;
    let d = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let t = {
        if vec3::dot(&d, &d) > 1.0e-20_f64.as_() {
            let ps = [
                p0[0] - point_pos[0],
                p0[1] - point_pos[1],
                p0[2] - point_pos[2],
            ];
            let a = vec3::dot(&d, &d);
            let b = vec3::dot(&d, &ps);
            let mut r: T = -b / a;
            r = if r < 0_f64.as_() { 0_f64.as_() } else { r };
            r = if r > 1_f64.as_() { 1_f64.as_() } else { r };
            r
        } else {
            0.5_f64.as_()
        }
    };
    [p0[0] + t * d[0], p0[1] + t * d[1], p0[2] + t * d[2]]
}
