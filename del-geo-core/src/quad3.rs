#[inline]
pub fn position<Real>(
    s0: Real,
    s1: Real,
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let w0 = (one - s0) * (one - s1);
    let w1 = s0 * (one - s1);
    let w2 = s0 * s1;
    let w3 = (one - s0) * s1;
    [
        w0 * q0[0] + w1 * q1[0] + w2 * q2[0] + w3 * q3[0],
        w0 * q0[1] + w1 * q1[1] + w2 * q2[1] + w3 * q3[1],
        w0 * q0[2] + w1 * q1[2] + w2 * q2[2] + w3 * q3[2],
    ]
}

#[inline]
pub fn ds0<Real>(
    _s0: Real,
    s1: Real,
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    [
        -(one - s1) * q0[0] + (one - s1) * q1[0] + s1 * q2[0] - s1 * q3[0],
        -(one - s1) * q0[1] + (one - s1) * q1[1] + s1 * q2[1] - s1 * q3[1],
        -(one - s1) * q0[2] + (one - s1) * q1[2] + s1 * q2[2] - s1 * q3[2],
    ]
}

#[inline]
pub fn ds1<Real>(
    s0: Real,
    _s1: Real,
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    [
        -(one - s0) * q0[0] - s0 * q1[0] + s0 * q2[0] + (one - s0) * q3[0],
        -(one - s0) * q0[1] - s0 * q1[1] + s0 * q2[1] + (one - s0) * q3[1],
        -(one - s0) * q0[2] - s0 * q1[2] + s0 * q2[2] + (one - s0) * q3[2],
    ]
}

#[inline]
pub fn ds01<Real>(
    _s0: Real,
    _s1: Real,
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> [Real; 3]
where
    Real: num_traits::Float,
{
    [
        q0[0] - q1[0] + q2[0] - q3[0],
        q0[1] - q1[1] + q2[1] - q3[1],
        q0[2] - q1[2] + q2[2] - q3[2],
    ]
}

pub fn nearest_to_origin<Real>(
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> ([Real; 3], Real, Real)
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let half = one / (one + one);
    let seeds = [
        (zero, zero),
        (one, zero),
        (one, one),
        (zero, one),
        (half, half),
    ];

    let mut dist_min = -one;
    let mut q_min = [zero; 3];
    let mut r0_min = zero;
    let mut r1_min = zero;

    for (mut r0, mut r1) in seeds {
        for _itr in 0..4 {
            let q = position(r0, r1, q0, q1, q2, q3);
            let dq0 = ds0(r0, r1, q0, q1, q2, q3);
            let dq1 = ds1(r0, r1, q0, q1, q2, q3);
            let ddq = ds01(r0, r1, q0, q1, q2, q3);

            let f0 = -crate::vec3::dot(&dq0, &q);
            let f1 = -crate::vec3::dot(&dq1, &q);

            let a00 = crate::vec3::dot(&dq0, &dq0);
            let a11 = crate::vec3::dot(&dq1, &dq1);
            let a01 = crate::vec3::dot(&dq1, &dq0) + crate::vec3::dot(&ddq, &q);

            let det = a00 * a11 - a01 * a01;
            let eps: Real = Real::from(1.0e-30f64).unwrap();
            if det.abs() < eps {
                break;
            }
            let inv_det = one / det;

            let dr0 = a11 * inv_det * f0 + (-a01 * inv_det) * f1;
            let dr1 = (-a01 * inv_det) * f0 + a00 * inv_det * f1;

            r0 = r0 + dr0;
            r1 = r1 + dr1;
        }

        let q = position(r0, r1, q0, q1, q2, q3);
        let tol: Real = Real::from(1.0e-4f64).unwrap();
        if r0 > -tol && r0 < one + tol && r1 > -tol && r1 < one + tol {
            let d = crate::vec3::norm(&q);
            if dist_min < zero || d < dist_min {
                dist_min = d;
                q_min = q;
                r0_min = r0;
                r1_min = r1;
            }
        }
    }

    if dist_min > zero {
        return (q_min, r0_min, r1_min);
    }

    {
        let (q, _, t) = crate::edge3::nearest_to_origin3(q0, q1);
        let d = crate::vec3::norm(&q);
        if dist_min < zero || d < dist_min {
            dist_min = d;
            q_min = q;
            r0_min = t;
            r1_min = zero;
        }
    }
    {
        let (q, _, t) = crate::edge3::nearest_to_origin3(q1, q2);
        let d = crate::vec3::norm(&q);
        if dist_min < zero || d < dist_min {
            dist_min = d;
            q_min = q;
            r0_min = one;
            r1_min = t;
        }
    }
    {
        let (q, _, t) = crate::edge3::nearest_to_origin3(q2, q3);
        let d = crate::vec3::norm(&q);
        if dist_min < zero || d < dist_min {
            dist_min = d;
            q_min = q;
            r0_min = one - t;
            r1_min = one;
        }
    }
    {
        let (q, _, t) = crate::edge3::nearest_to_origin3(q3, q0);
        let d = crate::vec3::norm(&q);
        if dist_min < zero || d < dist_min {
            let _ = dist_min;
            q_min = q;
            r0_min = zero;
            r1_min = one - t;
        }
    }

    (q_min, r0_min, r1_min)
}

pub fn intersection_against_line_bilinear<Real>(
    src: &[Real; 3],
    dir: &[Real; 3],
    q0: &[Real; 3],
    q1: &[Real; 3],
    q2: &[Real; 3],
    q3: &[Real; 3],
) -> Option<([Real; 3], Real, Real)>
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let half = one / (one + one);
    let (u, v) = crate::vec3::basis_xy_from_basis_z(dir);

    let seeds = [
        (zero, zero),
        (one, zero),
        (one, one),
        (zero, one),
        (half, half),
    ];

    let mut dist_min = -one;
    let mut q_min = [zero; 3];
    let mut r0_min = zero;
    let mut r1_min = zero;

    for (mut r0, mut r1) in seeds {
        for _itr in 0..4 {
            let q = position(r0, r1, q0, q1, q2, q3);
            let pq = crate::vec3::sub(&q, src);
            let dq0 = ds0(r0, r1, q0, q1, q2, q3);
            let dq1 = ds1(r0, r1, q0, q1, q2, q3);

            let f0 = -crate::vec3::dot(&u, &pq);
            let f1 = -crate::vec3::dot(&v, &pq);

            let a00 = crate::vec3::dot(&u, &dq0);
            let a01 = crate::vec3::dot(&u, &dq1);
            let a10 = crate::vec3::dot(&v, &dq0);
            let a11 = crate::vec3::dot(&v, &dq1);

            let det = a00 * a11 - a01 * a10;
            let eps: Real = Real::from(1.0e-30f64).unwrap();
            if det.abs() < eps {
                break;
            }
            let inv_det = one / det;

            let dr0 = a11 * inv_det * f0 + (-a01 * inv_det) * f1;
            let dr1 = (-a10 * inv_det) * f0 + a00 * inv_det * f1;

            r0 = r0 + dr0;
            r1 = r1 + dr1;
        }

        let q = position(r0, r1, q0, q1, q2, q3);
        let tol: Real = Real::from(1.0e-4f64).unwrap();
        if r0 > -tol && r0 < one + tol && r1 > -tol && r1 < one + tol {
            let d = crate::vec3::distance(&q, src);
            if dist_min < zero || d < dist_min {
                dist_min = d;
                q_min = q;
                r0_min = r0;
                r1_min = r1;
            }
        }
    }

    if dist_min > zero {
        Some((q_min, r0_min, r1_min))
    } else {
        None
    }
}
