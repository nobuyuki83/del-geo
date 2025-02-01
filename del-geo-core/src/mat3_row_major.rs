//! functions and the trait for 3x3 matrix with row major storage

/// trait for 3x3 matrix in row-major order
pub trait Mat3RowMajor<Real: num_traits::Float> {
    fn from_identity() -> Self;
    fn determinant(&self) -> Real;
    fn squared_norm(&self) -> Real;
    fn transpose(&self) -> Self;
    fn mult_mat_row_major(&self, b: &Self) -> Self;
    fn add(&self, b: &Self) -> Self;
    fn sub(&self, b: &Self) -> Self;
    fn scale(&self, s: Real) -> Self;
}
impl<Real> Mat3RowMajor<Real> for [Real; 9]
where
    Real: num_traits::Float + std::iter::Sum,
{
    fn from_identity() -> Self {
        from_identity()
    }
    fn determinant(&self) -> Real {
        determinant(self)
    }
    fn squared_norm(&self) -> Real {
        squared_norm(self)
    }
    fn transpose(&self) -> Self {
        transpose(self)
    }
    fn mult_mat_row_major(&self, b: &Self) -> Self {
        mult_mat_row_major(self, b)
    }
    fn add(&self, b: &Self) -> Self {
        add(self, b)
    }
    fn sub(&self, b: &Self) -> Self {
        sub(self, b)
    }
    fn scale(&self, s: Real) -> Self {
        scale(self, s)
    }
}

pub fn from_identity<T>() -> [T; 9]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, zero, one, zero, zero, zero, one]
}

pub fn from_diagonal<T>(d: &[T; 3]) -> [T; 9]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    [d[0], zero, zero, zero, d[1], zero, zero, zero, d[2]]
}

pub fn from_columns<T>(c0: &[T; 3], c1: &[T; 3], c2: &[T; 3]) -> [T; 9]
where
    T: num_traits::Float,
{
    [
        c0[0], c1[0], c2[0], c0[1], c1[1], c2[1], c0[2], c1[2], c2[2],
    ]
}

pub fn to_columns<T>(a: &[T; 9]) -> ([T; 3], [T; 3], [T; 3])
where
    T: num_traits::Float,
{
    ([a[0], a[3], a[6]], [a[1], a[4], a[7]], [a[2], a[5], a[8]])
}

// above: from method
// --------------------------

pub fn determinant<Real>(u: &[Real; 9]) -> Real
where
    Real: num_traits::Float,
{
    u[0] * u[4] * u[8] + u[3] * u[7] * u[2] + u[6] * u[1] * u[5]
        - u[0] * u[7] * u[5]
        - u[6] * u[4] * u[2]
        - u[3] * u[1] * u[8]
}

pub fn squared_norm<Real>(u: &[Real; 9]) -> Real
where
    Real: num_traits::Float + std::iter::Sum,
{
    u.iter().map(|&v| v * v).sum()
}

pub fn transpose<Real>(m: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

pub fn mult_mat_row_major<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    let mut r = [Real::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i * 3 + j] = r[i * 3 + j] + a[i * 3 + k] * b[k * 3 + j]
            }
        }
    }
    r
}

pub fn add<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    std::array::from_fn(|i| a[i] + b[i])
}

pub fn sub<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    std::array::from_fn(|i| a[i] - b[i])
}

pub fn scale<Real>(a: &[Real; 9], s: Real) -> [Real; 9]
where
    Real: num_traits::Float,
{
    std::array::from_fn(|i| a[i] * s)
}

/*
fn sort_eigen<Real>(g: &mut [Real; 3], v: &mut [Real; 9])
where
    Real: num_traits::Float,
{
    if g[1] > g[0] {
        g.swap(0, 1);
        v.swap(0, 1);
        v.swap(3, 4);
        v.swap(6, 7);
    }
    if g[2] > g[1] {
        g.swap(1, 2);
        v.swap(1, 2);
        v.swap(4, 5);
        v.swap(7, 8);
    }
    if g[1] > g[0] {
        g.swap(0, 1);
        v.swap(0, 1);
        v.swap(3, 4);
        v.swap(6, 7);
    }
}
 */

/*
/// m = UGV^T
pub fn svd_jacobi<Real>(m: &[Real; 9], nitr: usize) -> ([Real; 9], [Real; 9], [Real; 9])
where
    Real: num_traits::Float + std::iter::Sum + std::fmt::Debug,
{
    let one = Real::one();
    let zero = Real::zero();
    // M^TM = VGGV^T
    let mtm = [
        m[0] * m[0] + m[3] * m[3] + m[6] * m[6], // M^tM[0,0]
        m[1] * m[1] + m[4] * m[4] + m[7] * m[7], // M^tM[1,1]
        m[2] * m[2] + m[5] * m[5] + m[8] * m[8], // M^tM[2,2]
        m[1] * m[2] + m[4] * m[5] + m[7] * m[8], // M^tM[1,2]
        m[2] * m[0] + m[5] * m[3] + m[8] * m[6], // M^tM[2,0]
        m[0] * m[1] + m[3] * m[4] + m[6] * m[7], // M^tM[0,1]
    ];
    let Some((mut v, mut lv)) = crate::mat3_sym::eigen_decomposition_jacobi(&mtm, nitr) else {
        todo!()
    };
    assert!(lv[0] < lv[1] && lv[1] < lv[2]);
    sort_eigen(&mut lv, &mut v);
    dbg!(lv);
    // assert!(lv[0]<lv[1] && lv[1]<lv[2]);
    let mut g = lv.map(|x| x.max(zero).sqrt());

    let mut u0 = [
        m[0] * v[0] + m[1] * v[3] + m[2] * v[6],
        m[3] * v[0] + m[4] * v[3] + m[5] * v[6],
        m[6] * v[0] + m[7] * v[3] + m[8] * v[6],
    ];
    let mut u1 = [
        m[0] * v[1] + m[1] * v[4] + m[2] * v[7],
        m[3] * v[1] + m[4] * v[4] + m[5] * v[7],
        m[6] * v[1] + m[7] * v[4] + m[8] * v[7],
    ];
    let mut u2 = [
        m[0] * v[2] + m[1] * v[5] + m[2] * v[8],
        m[3] * v[2] + m[4] * v[5] + m[5] * v[8],
        m[6] * v[2] + m[7] * v[5] + m[8] * v[8],
    ];

    let mut u = [zero; 9];
    if v.determinant() < zero {
        // making right hand coordinate
        v[2] = -v[2]; // v[0,2] = v[0*3+2]
        v[5] = -v[5]; // v[1,2] = v[1*3+2]
        v[2 * 3 + 2] = -v[2 * 3 + 2];
        g[2] = -g[2];
        u[2] = -u[2]; // u[0,2] = u[0*3+2]
        u[5] = -u[5]; // u[1,2] = u[1*3+2]
        u[8] = -u[8]; // u[2,2] = u[2*3+2]
    }

    use crate::vec3::Vec3;
    let sql0 = u0.squared_norm();
    if sql0 > Real::epsilon() {
        u0.normalize_in_place();
        let sql1 = u1.squared_norm();
        if sql1 < Real::epsilon() {
            u1 = u0.map(|x| one - x.abs());
        } else {
            u1.normalize_in_place();
        }
        let d01 = u0[0] * u1[0] + u0[1] * u1[1] + u0[2] * u1[2];
        u1 = std::array::from_fn(|i| u1[i] - d01 * u0[i]);
        u1.normalize_in_place();
        let s2 = u0.cross(&u1);
        let d22 = u2[0] * s2[0] + u2[1] * s2[1] + u2[2] * s2[2];
        u2 = s2;
        if d22 < zero {
            g[2] = -g[2];
        }
    } else {
        u0 = [one, zero, zero];
        u1 = [zero, one, zero];
        u2 = [zero, zero, one];
    }
    u[0] = u0[0];
    u[1] = u1[0];
    u[2] = u2[0];
    u[3] = u0[1];
    u[4] = u1[1];
    u[5] = u2[1];
    u[6] = u0[2];
    u[7] = u1[2];
    u[8] = u2[2];
    let s = [g[0], zero, zero, zero, g[1], zero, zero, zero, g[2]];
    (u, s, v)
}
 */

/// Singluar Value Decomposition (SVD)
/// input = U * S * V^t
///
pub fn svd<Real>(
    f: &[Real; 9],
    mode: crate::mat3_sym::EigenDecompositionModes,
) -> Option<([Real; 9], [Real; 3], [Real; 9])>
where
    Real: num_traits::Float + std::iter::Sum + num_traits::FloatConst,
{
    let zero = Real::zero();
    let ft_f = f.transpose().mult_mat_row_major(f);
    let ft_f = crate::mat3_sym::from_mat3(&ft_f);
    let (v, mut lambda) = crate::mat3_sym::eigen_decomposition(&ft_f, mode)?;
    lambda.iter_mut().for_each(|x| *x = (*x).max(zero).sqrt());
    let ul = f.mult_mat_row_major(&v);
    let (mut u0, mut u1, mut u2) = to_columns(&ul);
    crate::vec3::normalize_in_place(&mut u0);
    crate::vec3::normalize_in_place(&mut u1);
    crate::vec3::normalize_in_place(&mut u2);
    let u = from_columns(&u0, &u1, &u2);
    Some((u, lambda, v))
}

pub fn enforce_rotation_matrix_for_svd<Real>(
    u: &[Real; 9],
    l: &[Real; 3],
    v: &[Real; 9],
) -> ([Real; 9], [Real; 3], [Real; 9])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    if determinant(v) < Real::zero() || determinant(u) < Real::zero() {
        let mut u = *u;
        let mut l = *l;
        let mut v = *v;
        if determinant(&v) < Real::zero() {
            v[2] = -v[2]; // v[0,2] = v[0*3+2]
            v[5] = -v[5]; // v[1,2] = v[1*3+2]
            v[8] = -v[8];
            l[2] = -l[2];
        }
        if determinant(&u) < Real::zero() {
            u[2] = -u[2]; // v[0,2] = v[0*3+2]
            u[5] = -u[5]; // v[1,2] = v[1*3+2]
            u[8] = -u[8];
            l[2] = -l[2];
        }
        (u, l, v)
    } else {
        (*u, *l, *v)
    }
}

#[test]
fn test_svd() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0);
    for (_iter, i_mode_eigen, is_rot) in itertools::iproduct!(0..100, 0..2, 0..2) {
        let m: [f64; 9] = std::array::from_fn(|_| rng.random_range(-1f64..1f64));
        let (u, s, v) = {
            let mode = match i_mode_eigen {
                0 => crate::mat3_sym::EigenDecompositionModes::JacobiNumIter(100),
                1 => crate::mat3_sym::EigenDecompositionModes::Analytic,
                _ => unreachable!(),
            };
            svd(&m, mode).unwrap()
        };
        let (u, s, v) = if is_rot == 1 {
            enforce_rotation_matrix_for_svd(&u, &s, &v)
        } else {
            (u, s, v)
        };
        if is_rot == 1 {
            let det_v = determinant(&v);
            assert!((det_v - 1.).abs() < 1.0e-10);
            let det_u = determinant(&u);
            assert!((det_u - 1.).abs() < 1.0e-10);
        }
        {
            // test u
            let diff = Mat3RowMajor::transpose(&u)
                .mult_mat_row_major(&u)
                .sub(&from_identity())
                .squared_norm();
            assert!(diff < 1.0e-20f64, "{}", diff);
        }
        {
            // test V V^t = I
            let diff = Mat3RowMajor::transpose(&v)
                .mult_mat_row_major(&v)
                .sub(&from_identity())
                .squared_norm();
            assert!(diff < 1.0e-20f64, "{}", diff);
        }
        {
            // test A = USV^t
            let s = from_diagonal(&s);
            let diff = u
                .mult_mat_row_major(&s)
                .mult_mat_row_major(&Mat3RowMajor::transpose(&v))
                .sub(&m)
                .squared_norm();
            assert!(diff < 1.0e-20f64, "{}", diff);
        }
    }
}

/// Jacobian of singular value decomposition
///
/// Papadopoulo, Théodore & Lourakis, Manolis. (2000). "
/// Estimating the Jacobian of the Singular Value Decomposition"
/// Theory and Applications. 554-570.
pub fn svd_differential(u: [f64; 9], s: [f64; 9], v: [f64; 9]) -> [[[f64; 3]; 3]; 9] {
    let inv_mat2 = |mut a0: f64, a1: f64| -> (f64, f64) {
        if (a0 - a1).abs() < 1.0e-6 {
            a0 += 1.0e-6;
        }
        let det_inv = 1.0 / (a0 * a0 - a1 * a1);
        (a0 * det_inv, -a1 * det_inv)
    };

    let ai0 = inv_mat2(s[4], s[8]);
    let ai1 = inv_mat2(s[8], s[0]);
    let ai2 = inv_mat2(s[0], s[4]);

    let mut diff = [[[0.0; 3]; 3]; 9];
    for (i, j) in itertools::iproduct!(0..3, 0..3) {
        // dSdu
        diff[3][i][j] = u[3 * i] * v[3 * j];
        diff[4][i][j] = u[3 * i + 1] * v[3 * j + 1];
        diff[5][i][j] = u[3 * i + 2] * v[3 * j + 2];
        {
            let b0 = [-u[3 * i + 2] * v[3 * j + 1], u[3 * i + 1] * v[3 * j + 2]];
            diff[0][i][j] = b0[0] * ai0.0 + b0[1] * ai0.1;
            diff[6][i][j] = -b0[0] * ai0.1 - b0[1] * ai0.0;
        }
        {
            let b1 = [-u[3 * i] * v[3 * j + 2], u[3 * i + 2] * v[3 * j]];
            diff[1][i][j] = b1[0] * ai1.0 + b1[1] * ai1.1;
            diff[7][i][j] = -b1[0] * ai1.1 - b1[1] * ai1.0;
        }
        {
            let b2 = [-u[3 * i + 1] * v[3 * j], u[3 * i] * v[3 * j + 1]];
            diff[2][i][j] = b2[0] * ai2.0 + b2[1] * ai2.1;
            diff[8][i][j] = -b2[0] * ai2.1 - b2[1] * ai2.0;
        }
    }
    diff
}

#[test]
fn test_svd_differential() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0);
    let eps = 1.0e-6;
    for _iter in 0..100 {
        let m0: [f64; 9] = std::array::from_fn(|_| rng.random::<f64>());
        let (u0, s0, v0) = svd(
            &m0,
            crate::mat3_sym::EigenDecompositionModes::JacobiNumIter(100),
        )
        .unwrap();
        let s0 = crate::mat3_row_major::from_diagonal(&s0);
        let diff = svd_differential(u0, s0, v0);
        for (i, j) in itertools::iproduct!(0..3, 0..3) {
            let m1 = {
                let mut m1 = m0;
                m1[i * 3 + j] += eps;
                m1
            };
            let (u1, s1, v1) = svd(
                &m1,
                crate::mat3_sym::EigenDecompositionModes::JacobiNumIter(100),
            )
            .unwrap();
            let s1 = from_diagonal(&s1);
            {
                let du = Mat3RowMajor::transpose(&u0)
                    .mult_mat_row_major(&u1.sub(&u0))
                    .scale(1. / eps);
                let v0 = diff[0][i][j];
                let v1 = diff[1][i][j];
                let v2 = diff[2][i][j];
                assert!(
                    (du[1 * 3 + 2] - v0).abs() < 1.0e-4 * (1.0 + v0.abs()),
                    "{}",
                    v0
                );
                assert!(
                    (du[2 * 3 + 0] - v1).abs() < 1.0e-4 * (1.0 + v1.abs()),
                    "{}",
                    v1
                );
                assert!(
                    (du[0 * 3 + 1] - v2).abs() < 1.0e-4 * (1.0 + v2.abs()),
                    "{}",
                    v2
                );
            }
            {
                let ds = [
                    (s1[0 * 3 + 0] - s0[0 * 3 + 0]) / eps,
                    (s1[1 * 3 + 1] - s0[1 * 3 + 1]) / eps,
                    (s1[2 * 3 + 2] - s0[2 * 3 + 2]) / eps,
                ];
                let v0 = diff[3][i][j];
                let v1 = diff[4][i][j];
                let v2 = diff[5][i][j];
                assert!((ds[0] - v0).abs() < 1.0e-5 * (1.0 + v0.abs()), "{}", v0);
                assert!((ds[1] - v1).abs() < 1.0e-5 * (1.0 + v1.abs()), "{}", v1);
                assert!((ds[2] - v2).abs() < 1.0e-5 * (1.0 + v2.abs()), "{}", v2);
            }
            {
                let dv = Mat3RowMajor::transpose(&v0)
                    .mult_mat_row_major(&v1.sub(&v0))
                    .scale(1. / eps);
                let v0 = diff[6][i][j];
                let v1 = diff[7][i][j];
                let v2 = diff[8][i][j];
                assert!(
                    (dv[1 * 3 + 2] - v0).abs() < 1.0e-4 * (1.0 + v0.abs()),
                    "{}",
                    v0
                );
                assert!(
                    (dv[2 * 3 + 0] - v1).abs() < 1.0e-4 * (1.0 + v1.abs()),
                    "{}",
                    v1
                );
                assert!(
                    (dv[0 * 3 + 1] - v2).abs() < 1.0e-4 * (1.0 + v2.abs()),
                    "{}",
                    v2
                );
            }
        }
    }
}
