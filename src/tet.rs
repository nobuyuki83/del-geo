//! methods for 3D tetrahedron

use num_traits::AsPrimitive;

pub fn volume_<T>(
    v1: &[T;3],
    v2: &[T;3],
    v3: &[T;3],
    v4: &[T;3]) -> T
where T: num_traits::Float + 'static + Copy,
    f64: num_traits::AsPrimitive<T>
{
    let a0 =  (v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]));
    let a1 = -(v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]));
    let a2 =  (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]));
    (a0 + a1 + a2) * 0.166_666_666_666_666_67_f64.as_()
}

pub fn volume<T>(
    v1: &nalgebra::Vector3<T>,
    v2: &nalgebra::Vector3<T>,
    v3: &nalgebra::Vector3<T>,
    v4: &nalgebra::Vector3<T>) -> T
    where T: nalgebra::RealField + 'static + Copy,
          f64: num_traits::AsPrimitive<T>
{
    let a0 =  (v2[0] - v1[0]) * ((v3[1] - v1[1]) * (v4[2] - v1[2]) - (v4[1] - v1[1]) * (v3[2] - v1[2]));
    let a1 = -(v2[1] - v1[1]) * ((v3[0] - v1[0]) * (v4[2] - v1[2]) - (v4[0] - v1[0]) * (v3[2] - v1[2]));
    let a2 =  (v2[2] - v1[2]) * ((v3[0] - v1[0]) * (v4[1] - v1[1]) - (v4[0] - v1[0]) * (v3[1] - v1[1]));
    (a0 + a1 + a2) * 0.166_666_666_666_666_67_f64.as_()
}


/**
 height of tetrahedron base:v1,v2,v3, top:v4
 */
pub fn height<T>(
    v1: &nalgebra::Vector3<T>,
    v2: &nalgebra::Vector3<T>,
    v3: &nalgebra::Vector3<T>,
    v4: &nalgebra::Vector3<T>) -> T
where T: nalgebra::RealField
{
    let na = (v2-v1).cross(&(v3-v1)).normalize();
    na.dot(&(v4-v1))
}