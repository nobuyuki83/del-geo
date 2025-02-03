//! methods for 3D tetrahedron

pub fn volume<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3], v4: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let three = T::one() + T::one() + T::one();
    let one_6th = T::one() / (three + three);
    let a0 =
        (v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]));
    let a1 =
        -(v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]));
    let a2 =
        (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]));
    (a0 + a1 + a2) * one_6th
}
