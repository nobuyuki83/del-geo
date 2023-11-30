//! methods for 3D plane

use num_traits::AsPrimitive;

pub fn intersection_line3<T>(
    o: &nalgebra::Vector3<T>, // one point on plane
    n: &nalgebra::Vector3<T>, // plane normal
    s: &nalgebra::Vector3<T>, // one point on line
    d: &nalgebra::Vector3<T>) // direction of line
    -> nalgebra::Vector3<T>
    where T: nalgebra::RealField
{
    let t = ((o - s).dot(n)) / (d.dot(n));
    s + d.scale(t)
}

pub fn intersection_ray3<T>(
    o: &nalgebra::Vector3<T>, // one point on plane
    n: &nalgebra::Vector3<T>, // plane normal
    s: &nalgebra::Vector3<T>, // one point on line
    d: &nalgebra::Vector3<T>) // direction of line
    -> Option<T>
    where T: nalgebra::RealField
{
    let t = ((o - s).dot(n)) / (d.dot(n));
    if t < T::zero() { None } else { Some(t)}
}


pub fn nearest_to_point3<T>(
    p: &nalgebra::Vector3<T>, // point
    o: &nalgebra::Vector3<T>, // origin
    n: &nalgebra::Vector3<T>) -> nalgebra::Vector3<T> // normal
    where T: nalgebra::RealField
{
    let n0 = n.normalize();
    p + n0.scale((o - p).dot(&n0))
}

pub fn intersection_line3_triplane3<T>(
    src: &nalgebra::Vector3<T>,
    dir: &nalgebra::Vector3<T>,
    q0: &nalgebra::Vector3<T>,
    q1: &nalgebra::Vector3<T>,
    q2: &nalgebra::Vector3<T>,
    eps: T) -> Option<(nalgebra::Vector3<T>, T, T, T)>
where T: nalgebra::RealField + 'static + Copy,
      f64: num_traits::AsPrimitive<T>
{
    let mut r0 = crate::tet::volume(src, &(src + dir), q1, q2);
    let mut r1 = crate::tet::volume(src, &(src + dir), q2, q0);
    let mut r2 = crate::tet::volume(src, &(src + dir), q0, q1);
    let v012 = r0 + r1 + r2;
    let v012_inv = 1_f64.as_() / v012;
    r0 *= v012_inv;
    r1 *= v012_inv;
    r2 *= v012_inv;
    let p0 = q0.scale(r0) + q1.scale(r1) + q2.scale(r2);
    if r0 > eps && r1 > eps && r2 > eps {
        return Some((p0, r0, r1, r2));
    }
    None
}