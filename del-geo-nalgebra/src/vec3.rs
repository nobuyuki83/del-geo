use num_traits::AsPrimitive;

pub fn scalar_triple_product<T>(
    a: &nalgebra::Vector3<T>,
    b: &nalgebra::Vector3<T>,
    c: &nalgebra::Vector3<T>,
) -> T
where
    T: nalgebra::RealField,
{
    b.cross(c).dot(a)
}

pub fn frame_from_z_vector<T>(
    vec_n: nalgebra::Vector3<T>,
) -> (nalgebra::Vector3<T>, nalgebra::Vector3<T>)
where
    T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>,
{
    let vec_s = nalgebra::Vector3::<T>::new(T::zero(), T::one(), T::zero());
    let mut vec_x = vec_s.cross(&vec_n);
    let len = vec_x.norm();
    if len < 1.0e-10_f64.as_() {
        let vec_t = nalgebra::Vector3::<T>::new(T::one(), T::zero(), T::zero());
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

pub fn sample_unit_cube<T>() -> nalgebra::Vector3<T>
where
    T: nalgebra::RealField + nalgebra::Scalar,
    rand::distributions::Standard: rand::prelude::Distribution<T>,
{
    use rand::Rng;
    let mut p0 = nalgebra::Vector3::<T>::zeros();
    let mut rng = rand::thread_rng();
    for v in p0.iter_mut() {
        *v = rng.gen();
    }
    p0
}

pub fn from_homogeneous<T>(v: &nalgebra::Vector4<T>) -> Option<nalgebra::Vector3<T>>
where
    T: Copy + nalgebra::RealField,
{
    if v[3].is_zero() {
        return None;
    }
    Some(nalgebra::Vector3::<T>::new(
        v[0] / v[3],
        v[1] / v[3],
        v[2] / v[3],
    ))
}
