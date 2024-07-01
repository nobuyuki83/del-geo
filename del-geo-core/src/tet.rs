//! methods for 3D tetrahedron

use num_traits::AsPrimitive;

pub fn volume_<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3], v4: &[T; 3]) -> T
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let a0 =
        (v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]));
    let a1 =
        -(v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]));
    let a2 =
        (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]));
    (a0 + a1 + a2) * 0.166_666_666_666_666_67_f64.as_()
}
