//! methods for 3D vector

use num_traits::AsPrimitive;

pub fn squared_norm_<T>(p: &[T]) -> T
    where T: std::ops::Mul<Output=T> + std::ops::Add<Output=T> + Copy
{
    assert_eq!(p.len(), 3);
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
}

pub fn norm_<T>(
    v: &[T]) -> T
    where T: num_traits::Float
{
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn normalize_<T>(
    v: &mut [T]) -> T
where T: num_traits::Float + std::ops::MulAssign
{
    let l = norm_(v);
    let linv = T::one() / l;
    v[0] *= linv;
    v[1] *= linv;
    v[2] *= linv;
    l
}

pub fn cross_mut_<T>(
    vo: &mut [T],
    v1: &[T],
    v2: &[T])
    where T: std::ops::Mul<Output=T> + std::ops::Sub<Output=T> + Copy
{
    vo[0] = v1[1] * v2[2] - v2[1] * v1[2];
    vo[1] = v1[2] * v2[0] - v2[2] * v1[0];
    vo[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

pub fn cross_<T>(
    v1: &[T],
    v2: &[T]) -> [T; 3]
    where T: std::ops::Mul<Output=T> + std::ops::Sub<Output=T> + Copy
{
    [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]]
}

pub fn dot_<T>(
    a: &[T],
    b: &[T]) -> T
    where T: std::ops::Mul<Output=T> + std::ops::Add<Output=T> + Copy
{
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn sub_<T>(
    a: &[T],
    b: &[T]) -> [T; 3]
    where T: std::ops::Sub<Output=T> + Copy
{
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale_<T>(
    a: &[T],
    s: T) -> [T; 3]
    where T: Copy + std::ops::Mul<Output=T> {
    [s * a[0], s * a[1], s * a[2]]
}

pub fn distance_<T>(
    p0: &[T],
    p1: &[T]) -> T
    where T: num_traits::Float
{
    let v0 = p1[0] - p0[0];
    let v1 = p1[1] - p0[1];
    let v2 = p1[2] - p0[2];
    (v0 * v0 + v1 * v1 + v2 * v2).sqrt()
}

pub fn scalar_triple_product_<T>(
    a: &[T],
    b: &[T],
    c: &[T]) -> T
    where T: std::ops::Mul<Output=T> + std::ops::Sub<Output=T> + std::ops::Add<Output=T> + Copy
{
    let v0: T = a[0] * (b[1] * c[2] - b[2] * c[1]);
    let v1: T = a[1] * (b[2] * c[0] - b[0] * c[2]);
    let v2: T = a[2] * (b[0] * c[1] - b[1] * c[0]);
    v0 + v1 + v2
}

// -------------------

pub fn scalar_triple_product<T>(
    a: &nalgebra::Vector3::<T>,
    b: &nalgebra::Vector3::<T>,
    c: &nalgebra::Vector3::<T>) -> T
    where T: nalgebra::RealField
{
    b.cross(c).dot(a)
}

pub fn frame_from_z_vector<T>(
    vec_n: nalgebra::Vector3::<T>) -> (nalgebra::Vector3::<T>, nalgebra::Vector3::<T>)
where T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>
{
    let vec_s = nalgebra::Vector3::<T>::new(T::zero(), T::one(), T::zero() );
    let mut vec_x = vec_s.cross(&vec_n);
    let len = vec_x.norm();
    if len < 1.0e-10_f64.as_() {
        let vec_t = nalgebra::Vector3::<T>::new( T::one(), T::zero(), T::zero() );
        let vec_x = vec_t.cross(&vec_n);
        let vec_y = vec_n.cross(&vec_x);
        (vec_x, vec_y)
    } else {
        let invlen = T::one() / len;
        vec_x *= invlen;
        let vec_y = vec_n.cross(&vec_x);
        (vec_x, vec_y)
    }
}

pub fn sample_unit_cube<T>() -> nalgebra::Vector3::<T>
where T: nalgebra::RealField + nalgebra::Scalar,
      rand::distributions::Standard: rand::prelude::Distribution<T>
{
    use rand::Rng;
    let mut p0 = nalgebra::Vector3::<T>::zeros();
    let mut rng = rand::thread_rng();
    for v in p0.iter_mut() {
        *v = rng.gen();
    }
    p0
}

pub fn to_na<T>(vtx2xyz: &[T], i_vtx: usize) -> nalgebra::Vector3::<T>
    where T: Copy + nalgebra::RealField
{
    nalgebra::Vector3::<T>::from_row_slice(&vtx2xyz[i_vtx *3..(i_vtx +1)*3])
}