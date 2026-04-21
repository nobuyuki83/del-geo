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

fn barycentric_coord_for_origin<Real>(
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> Option<(Real, Real, Real)>
where
    Real: num_traits::Float,
{
    let total = volume(q0, q1, q2, q3);
    let zero = Real::zero();
    let eps = Real::epsilon();
    if total.abs() < eps {
        return None;
    }
    let origin = &[zero; 3];
    let r0 = volume(origin, q1, q2, q3) / total;
    let r1 = volume(q0, origin, q2, q3) / total;
    let r2 = volume(q0, q1, origin, q3) / total;
    Some((r0, r1, r2))
}

#[allow(unused_assignments)]
pub fn nearest_to_origin<Real>(
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> ([Real; 3], Real, Real, Real)
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let zero = Real::zero();
    let one = Real::one();
    let mut r0 = zero;
    let mut r1 = zero;
    let mut r2 = zero;
    let mut r3 = zero;

    let mut p_min = *q0;

    // inside tetrahedron
    if let Some((rr0, rr1, rr2)) = barycentric_coord_for_origin(q0, q1, q2, q3) {
        r0 = rr0;
        r1 = rr1;
        r2 = rr2;
        r3 = one - r0 - r1 - r2;
        p_min = crate::vec3::add4(q0, r0, q1, r1, q2, r2, q3, r3);
        if r0 > zero && r1 > zero && r2 > zero && r3 > zero {
            return (p_min, r0, r1, r2);
        }
    }

    // faces
    {
        // face123
        let (p, a1, a2, a3) = crate::tri3::nearest_to_origin3(q1, q2, q3);
        p_min = p;
        r0 = zero;
        r1 = a1;
        r2 = a2;
        r3 = a3;
        if r1 > zero && r2 > zero && r3 > zero {
            return (p_min, r0, r1, r2);
        }
    }
    {
        // face230
        let (p, a2, a3, a0) = crate::tri3::nearest_to_origin3(q2, q3, q0);
        p_min = p;
        r0 = a0;
        r1 = zero;
        r2 = a2;
        r3 = a3;
        if r2 > zero && r3 > zero && r0 > zero {
            return (p_min, r0, r1, r2);
        }
    }
    {
        // face301
        let (p, a3, a0, a1) = crate::tri3::nearest_to_origin3(q3, q0, q1);
        p_min = p;
        r0 = a0;
        r1 = a1;
        r2 = zero;
        r3 = a3;
        if r3 > zero && r0 > zero && r1 > zero {
            return (p_min, r0, r1, r2);
        }
    }
    {
        // face012
        let (p, a0, a1, a2) = crate::tri3::nearest_to_origin3(q0, q1, q2);
        p_min = p;
        r0 = a0;
        r1 = a1;
        r2 = a2;
        r3 = zero;
        if r0 > zero && r1 > zero && r2 > zero {
            return (p_min, r0, r1, r2);
        }
    }

    // edges
    let mut d_min = crate::vec3::norm(q0);

    {
        // edge01
        let (p01, s0, s1) = crate::edge3::nearest_to_origin3(q0, q1);
        let d01 = crate::vec3::norm(&p01);
        if d01 < d_min {
            d_min = d01;
            p_min = p01;
            r0 = s0;
            r1 = s1;
            r2 = zero;
            r3 = zero;
        }
    }
    {
        // edge02
        let (p02, s0, s2) = crate::edge3::nearest_to_origin3(q0, q2);
        let d02 = crate::vec3::norm(&p02);
        if d02 < d_min {
            d_min = d02;
            p_min = p02;
            r0 = s0;
            r1 = zero;
            r2 = s2;
            r3 = zero;
        }
    }
    {
        // edge03
        let (p03, s0, s3) = crate::edge3::nearest_to_origin3(q0, q3);
        let d03 = crate::vec3::norm(&p03);
        if d03 < d_min {
            d_min = d03;
            p_min = p03;
            r0 = s0;
            r1 = zero;
            r2 = zero;
            r3 = s3;
        }
    }
    {
        // edge12
        let (p12, s1, s2) = crate::edge3::nearest_to_origin3(q1, q2);
        let d12 = crate::vec3::norm(&p12);
        if d12 < d_min {
            d_min = d12;
            p_min = p12;
            r0 = zero;
            r1 = s1;
            r2 = s2;
            r3 = zero;
        }
    }
    {
        // edge13
        let (p13, s1, s3) = crate::edge3::nearest_to_origin3(q1, q3);
        let d13 = crate::vec3::norm(&p13);
        if d13 < d_min {
            d_min = d13;
            p_min = p13;
            r0 = zero;
            r1 = s1;
            r2 = zero;
            r3 = s3;
        }
    }
    {
        // edge23
        let (p23, s2, s3) = crate::edge3::nearest_to_origin3(q2, q3);
        let d23 = crate::vec3::norm(&p23);
        if d23 < d_min {
            let _ = d_min;
            p_min = p23;
            r0 = zero;
            r1 = zero;
            r2 = s2;
            r3 = s3;
        }
    }

    (p_min, r0, r1, r2)
}

#[test]
fn test_nearest_to_origin() {
    // regular tet centered at origin — origin is strictly inside
    let q0 = [1.0f64, 1.0, -1.0];
    let q1 = [1.0, -1.0, 1.0];
    let q2 = [-1.0, 1.0, 1.0];
    let q3 = [-1.0, -1.0, -1.0];

    {
        let (p, r0, r1, r2) = nearest_to_origin(&q0, &q1, &q2, &q3);
        let r3 = 1.0 - r0 - r1 - r2;
        assert!(
            r0 > 0.0 && r1 > 0.0 && r2 > 0.0 && r3 > 0.0,
            "bary coords should all be positive inside"
        );
        assert!(
            crate::vec3::norm(&p) < 1.0e-10,
            "nearest point should be origin"
        );
        let recon = crate::vec3::add4(&q0, r0, &q1, r1, &q2, r2, &q3, r3);
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction mismatch"
        );
    }

    // tet shifted far from origin — origin is outside, nearest point is on a face
    {
        let q0 = [3.0f64, 0.0, 0.0];
        let q1 = [4.0, 1.0, 0.0];
        let q2 = [4.0, 0.0, 1.0];
        let q3 = [4.0, 0.0, 0.0];
        let (p, r0, r1, r2) = nearest_to_origin(&q0, &q1, &q2, &q3);
        let r3 = 1.0 - r0 - r1 - r2;
        assert!(
            (r0 + r1 + r2 + r3 - 1.0).abs() < 1.0e-10,
            "bary coords should sum to 1"
        );
        assert!(
            r0 >= -1.0e-10 && r1 >= -1.0e-10 && r2 >= -1.0e-10 && r3 >= -1.0e-10,
            "bary coords non-negative"
        );
        let d = crate::vec3::norm(&p);
        assert!(
            d <= crate::vec3::norm(&q0) + 1.0e-10,
            "nearest point not closer than q0"
        );
        assert!(
            d <= crate::vec3::norm(&q1) + 1.0e-10,
            "nearest point not closer than q1"
        );
        assert!(
            d <= crate::vec3::norm(&q2) + 1.0e-10,
            "nearest point not closer than q2"
        );
        assert!(
            d <= crate::vec3::norm(&q3) + 1.0e-10,
            "nearest point not closer than q3"
        );
        let recon = crate::vec3::add4(&q0, r0, &q1, r1, &q2, r2, &q3, r3);
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction mismatch"
        );
    }
}
