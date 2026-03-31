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

/// Ω(r1,r2,r3): 原点から見た三角形(r1,r2,r3)の有向立体角
pub fn solid_angle<T>(r1: &[T; 3], r2: &[T; 3], r3: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let zero = T::zero();
    let two = T::one() + T::one();
    let n1 = r1.norm();
    let n2 = r2.norm();
    let n3 = r3.norm();

    // 退化（原点と一致）を避ける：用途に応じてResultにしてもOK
    if n1 == zero || n2 == zero || n3 == zero {
        return zero;
    }

    let triple = r1.dot(&crate::vec3::cross(r2, r3)); // det[r1 r2 r3]
    let denom = n1 * n2 * n3 + r1.dot(r2) * n3 + r2.dot(r3) * n1 + r3.dot(r1) * n2;

    two * triple.atan2(denom)
}

pub fn gauss_linking_number_edge_edge<T>(a: &[T; 3], b: &[T; 3], c: &[T; 3], d: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    // Φ = Ω(c-a, d-a, b-a) - Ω(c-b, d-b, a-b)
    use crate::vec3::Vec3;
    let t1 = solid_angle(&c.sub(a), &d.sub(a), &b.sub(a));
    let t2 = solid_angle(&c.sub(b), &d.sub(b), &a.sub(b));
    t1 - t2
}
