//! methods for 3D plane

pub fn intersection_line3<T>(
    o: &[T; 3], // one point on plane
    n: &[T; 3], // plane normal
    s: &[T; 3], // one point on line
    d: &[T; 3],
) -> [T; 3]
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let t = o.sub(s).dot(n) / d.dot(n);
    s.add(&d.scale(t))
}

pub fn intersection_ray3<T>(
    o: &[T; 3], // one point on a plane
    n: &[T; 3], // plane normal
    s: &[T; 3], // one point on a line
    d: &[T; 3],
) -> Option<T>
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let t = o.sub(s).dot(n) / d.dot(n);
    if t < T::zero() {
        None
    } else {
        Some(t)
    }
}

pub fn nearest_to_point3<T>(
    p: &[T; 3], // point
    o: &[T; 3], // origin
    n: &[T; 3],
) -> [T; 3]
// normal
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let n0 = n.normalize();
    p.add(&n0.scale(o.sub(p).dot(&n0)))
}
