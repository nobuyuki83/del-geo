//! Symmetric 3x3 matrix class
//! the 3x3 matrix is stored in a 6-dim array
//! `[m[0,0], m[1,1], m[2,2], m[1,2], m[2,0], m[0,1]]`
//! (first diagonal, then off-diagonal)

/// the input can be both column major or row major
pub fn from_mat3_by_symmetrization<Real>(m: &[Real; 9]) -> [Real; 6]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let half = one / (one + one);
    [
        m[0],
        m[4],
        m[8],
        (m[5] + m[7]) * half,
        (m[2] + m[6]) * half,
        (m[1] + m[3]) * half,
    ]
}

pub fn to_mat3_row_major<Real>(sm: &[Real; 6]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    [
        sm[0], sm[5], sm[4], sm[5], sm[1], sm[3], sm[4], sm[3], sm[2],
    ]
}

/// Frobenius norm squared
pub fn squared_norm<Real>(sm: &[Real; 6]) -> Real
where
    Real: num_traits::Float,
{
    let two = Real::one() + Real::one();
    sm[0] * sm[0]
        + sm[1] * sm[1]
        + sm[2] * sm[2]
        + two * (sm[3] * sm[3] + sm[4] * sm[4] + sm[5] * sm[5])
}

pub fn trace<Real>(sm: &[Real; 6]) -> Real
where
    Real: num_traits::Float,
{
    sm[0] + sm[1] + sm[2]
}

pub fn determinant<Real>(sm: &[Real; 6]) -> Real
where
    Real: num_traits::Float,
{
    let two = Real::one() + Real::one();
    sm[0] * sm[1] * sm[2] + two * sm[3] * sm[4] * sm[5]
        - sm[0] * sm[3] * sm[3]
        - sm[1] * sm[4] * sm[4]
        - sm[2] * sm[5] * sm[5]
}

/// this function returns U and \Sigma
/// A = U * \Sigma * U^t
pub fn eigen_decomposition_jacobi<Real>(
    sm: &[Real; 6],
    num_iter: usize,
) -> Option<([Real; 9], [Real; 3])>
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    let half = one / two;
    let mut u = [zero; 9];
    // initialize u as identity matrix
    u[0] = one;
    u[4] = one;
    u[8] = one;
    let dnrm = squared_norm(sm);
    if dnrm < Real::epsilon() {
        return None;
    } // this matrix is too small
    let scale = dnrm.sqrt();
    let inv_scale = one / scale;
    let mut sms = sm.map(|x| x * inv_scale);

    for _itr in 0..num_iter {
        let m = sms;
        let v = u;
        let a12 = sms[3].abs();
        let a20 = sms[4].abs();
        let a01 = sms[5].abs();
        if a12 >= a20 && a12 >= a01 {
            // when a12 sms[3] is the biggest
            let t = half * (two * m[3]).atan2(m[2] - m[1]);
            let ct = t.cos();
            let st = t.sin();
            sms[1] = ct * ct * m[1] + st * st * m[2] - two * st * ct * m[3];
            sms[2] = ct * ct * m[2] + st * st * m[1] + two * st * ct * m[3];
            sms[3] = zero; // (ct*ct-st*st)*m[3]+st*ct*(m[1]-m[2]);
            sms[4] = st * m[5] + ct * m[4];
            sms[5] = ct * m[5] - st * m[4];
            //
            u[1] = ct * v[1] - st * v[2];
            u[2] = st * v[1] + ct * v[2];
            u[4] = ct * v[4] - st * v[5];
            u[5] = st * v[4] + ct * v[5];
            u[7] = ct * v[7] - st * v[8];
            u[8] = st * v[7] + ct * v[8];
        } else if a20 >= a01 && a20 >= a12 {
            // when a20 sms[4] is the biggest
            // the above condition statement shoud pass exactly once for each iteration.
            let t = half * (two * m[4]).atan2(m[2] - m[0]);
            let ct = t.cos();
            let st = t.sin();
            sms[0] = ct * ct * m[0] + st * st * m[2] - two * st * ct * m[4];
            sms[2] = ct * ct * m[2] + st * st * m[0] + two * st * ct * m[4];
            sms[3] = st * m[5] + ct * m[3];
            sms[4] = zero; // (ct*ct-st*st)*m[4]+st*ct*(m[0]-m[2]);
            sms[5] = ct * m[5] - st * m[3];
            //
            u[0] = ct * v[0] - st * v[2];
            u[2] = st * v[0] + ct * v[2];
            u[3] = ct * v[3] - st * v[5];
            u[5] = st * v[3] + ct * v[5];
            u[6] = ct * v[6] - st * v[8];
            u[8] = st * v[6] + ct * v[8];
        } else {
            // when a01 sms[5] is the biggest
            // the condition statement shoud pass exactly once for each iteration.
            let t = half * (two * m[5]).atan2(m[1] - m[0]);
            let ct = t.cos();
            let st = t.sin();
            sms[0] = ct * ct * m[0] + st * st * m[1] - two * st * ct * m[5];
            sms[1] = ct * ct * m[1] + st * st * m[0] + two * st * ct * m[5];
            sms[3] = st * m[4] + ct * m[3];
            sms[4] = ct * m[4] - st * m[3];
            sms[5] = zero; // (ct*ct-st*st)*m[5]+st*ct*(m[0]-m[1]);
            //
            u[0] = ct * v[0] - st * v[1];
            u[1] = st * v[0] + ct * v[1];
            u[3] = ct * v[3] - st * v[4];
            u[4] = st * v[3] + ct * v[4];
            u[6] = ct * v[6] - st * v[7];
            u[7] = st * v[6] + ct * v[7];
        }
    }
    let l = std::array::from_fn(|i| scale * sms[i]);
    Some((u, l))
}

pub fn eigen_values_analytic<Real>(m: &[Real; 6]) -> Option<[Real; 3]>
where
    Real: num_traits::Float + num_traits::FloatConst,
{
    let one = Real::one();
    let two = one + one;
    let three = two + one;
    let half = one / two;
    let one_third = one / three;
    let one_sixth = one_third * half;
    let p1 = m[5].powi(2) + m[4].powi(2) + m[3].powi(2);
    let q = trace(m) * one_third;
    let p2 = (m[0] - q).powi(2) + (m[1] - q).powi(2) + (m[2] - q).powi(2) + two * p1;
    let p = (p2 * one_sixth).sqrt();
    if p.abs() < Real::epsilon() {
        None
    } else {
        let b = {
            let s = one / p;
            [
                s * (m[0] - q),
                s * (m[1] - q),
                s * (m[2] - q),
                s * m[3],
                s * m[4],
                s * m[5],
            ]
        };
        let r = determinant(&b) * half;
        let phi = if r <= -one {
            Real::PI() * one_third
        } else if r >= one {
            Real::zero()
        } else {
            r.acos() * one_third
        };
        let eig_large = q + two * p * phi.cos();
        let eig_small = q + two * p * (phi + two * Real::PI() * one_third).cos();
        let eig_medium = three * q - eig_large - eig_small;
        Some([eig_small, eig_medium, eig_large])
    }
}

fn vector_with_largest_norm<Real>(a: &[Real; 3], b: &[Real; 3], c: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    use crate::vec3::Vec3;
    let (a_norm, b_norm, c_norm) = (a.squared_norm(), b.squared_norm(), c.squared_norm());
    if a_norm > b_norm {
        if a_norm > c_norm { *a } else { *c }
    } else if b_norm > c_norm {
        *b
    } else {
        *c
    }
}
fn find_orthogonal<Real>(sm: &[Real; 6]) -> ([Real; 3], Option<[Real; 3]>)
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let (u, v, w) = (
        [sm[0], sm[5], sm[4]],
        [sm[5], sm[1], sm[3]],
        [sm[4], sm[3], sm[2]],
    );
    use crate::vec3::Vec3;
    let (uv, vw, wu) = (u.cross(&v), v.cross(&w), w.cross(&u));
    let q = vector_with_largest_norm(&uv, &vw, &wu);
    if q.squared_norm() < Real::epsilon() {
        let p = vector_with_largest_norm(&u, &v, &w);
        let x = p.cross(&[one, zero, zero]);
        let x = if x.squared_norm() < Real::epsilon() {
            p.cross(&[zero, one, zero])
        } else {
            x
        };
        let y = p.cross(&x);
        (x.normalize(), Some(y.normalize()))
    } else {
        (q.normalize(), None)
    }
}

#[allow(non_snake_case)]
pub fn eigen_vectors_from_eigen_values<Real>(A: &[Real; 6], lmd: &[Real; 3]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    use crate::vec3::Vec3;
    let tmp0 = [
        A[0] - lmd[0],
        A[1] - lmd[0],
        A[2] - lmd[0],
        A[3],
        A[4],
        A[5],
    ];
    let uv = find_orthogonal(&tmp0);
    if let Some(v) = uv.1 {
        let c0 = uv.0;
        let c1 = v;
        let c2 = c0.cross(&c1);
        crate::mat3_row_major::from_columns(&c0, &c1, &c2)
    } else {
        let tmp1 = [
            A[0] - lmd[1],
            A[1] - lmd[1],
            A[2] - lmd[1],
            A[3],
            A[4],
            A[5],
        ];
        let tmp = find_orthogonal(&tmp1);
        let c0 = uv.0;
        let c1 = tmp.0;
        let c2 = c0.cross(&c1);
        crate::mat3_row_major::from_columns(&c0, &c1, &c2)
    }
}

/// this function returns U and \Sigma
/// A = U * \Sigma * U^t
pub fn eigen_decomposition_analytic<Real>(sm: &[Real; 6]) -> Option<([Real; 9], [Real; 3])>
where
    Real: num_traits::Float + num_traits::FloatConst,
{
    let l = eigen_values_analytic(sm)?;
    let u = eigen_vectors_from_eigen_values(sm, &l);
    Some((u, l))
}

pub enum EigenDecompositionModes {
    Analytic,
    JacobiNumIter(usize),
}

pub fn eigen_decomposition<Real>(
    sm: &[Real; 6],
    mode: EigenDecompositionModes,
) -> Option<([Real; 9], [Real; 3])>
where
    Real: num_traits::Float + num_traits::FloatConst,
{
    match mode {
        EigenDecompositionModes::Analytic => eigen_decomposition_analytic(sm),
        EigenDecompositionModes::JacobiNumIter(num_iter) => {
            eigen_decomposition_jacobi(sm, num_iter)
        }
    }
}

/*
       {
           let Some((_u, l_num)) = eigen_decomposition_jacobi(&sm, 20) else {
               todo!()
           };
           assert!((l_num[0] - l[0]).abs() < 1.0e-12);
           assert!((l_num[1] - l[1]).abs() < 1.0e-12);
           assert!((l_num[2] - l[2]).abs() < 1.0e-12);
       }
*/

#[test]
fn test_eigen_decomposition() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    // std::uniform_real_distribution < double > dist(-50.0, 50.0);
    for (_itr, i_mode) in itertools::iproduct!(0..10000, 0..2) {
        let sm = std::array::from_fn(|_| rng.random_range(-30f64..30f64));
        let (u, l) = match i_mode {
            0 => eigen_decomposition(&sm, EigenDecompositionModes::JacobiNumIter(100)).unwrap(),
            1 => eigen_decomposition(&sm, EigenDecompositionModes::Analytic).unwrap(),
            _ => unreachable!(),
        };
        assert!(l[0] < l[1] && l[1] < l[2]);
        use crate::mat3_row_major::Mat3RowMajor;
        let ut = u.transpose();
        {
            let utu = u.mult_mat_row_major(&ut);
            let err = crate::mat3_row_major::from_identity()
                .sub(&utu)
                .squared_norm();
            assert!(err < 1.0e-20);
        }
        {
            // check determinant of `u` is one
            let det = u.determinant();
            assert!((det - 1.).abs() < 1.0e-10, "{}", (det - 1.0).abs());
        }
        let l = crate::mat3_row_major::from_diagonal(&l);
        let ulut = u.mult_mat_row_major(&l).mult_mat_row_major(&ut);
        let sm = to_mat3_row_major(&sm);
        let err = sm.sub(&ulut).squared_norm();
        assert!(err < 1.0e-10);
    }
}
