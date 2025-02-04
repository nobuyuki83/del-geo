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

pub fn intersection_line3_triplane3<T>(
    src: &[T; 3],
    dir: &[T; 3],
    q0: &[T; 3],
    q1: &[T; 3],
    q2: &[T; 3],
    eps: T,
) -> Option<([T; 3], T, T, T)>
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let mut r0 = crate::tet::volume(src, &src.add(dir), q1, q2);
    let mut r1 = crate::tet::volume(src, &src.add(dir), q2, q0);
    let mut r2 = crate::tet::volume(src, &src.add(dir), q0, q1);
    let v012 = r0 + r1 + r2;
    let v012_inv = T::one() / v012;
    r0 = r0 * v012_inv;
    r1 = r1 * v012_inv;
    r2 = r2 * v012_inv;
    let p0 = crate::vec3::add_three(&q0.scale(r0), &q1.scale(r1), &q2.scale(r2));
    if r0 > eps && r1 > eps && r2 > eps {
        return Some((p0, r0, r1, r2));
    }
    None
}
