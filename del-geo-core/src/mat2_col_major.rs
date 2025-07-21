//! methods for 2x2 matrix (column major storage)

pub trait Mat2ColMajor<T: num_traits::Float>
where
    Self: Sized,
{
    fn scale(&self, s: T) -> Self;
    fn mult_mat_col_major(&self, other: &Self) -> Self;
    fn mult_vec(&self, other: &[T; 2]) -> [T; 2];
    fn add(&self, other: &Self) -> Self;
    fn add_in_place(&mut self, other: &Self);
    fn sub(&self, other: &Self) -> Self;
    fn transpose(&self) -> Self;
    fn squared_norm(&self) -> T;
    fn determinant(&self) -> T;
    fn set_zero(&mut self);
}

impl<Real> Mat2ColMajor<Real> for [Real; 4]
where
    Real: num_traits::Float,
{
    fn scale(&self, s: Real) -> Self {
        std::array::from_fn(|i| self[i] * s)
    }

    fn add(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i] + other[i])
    }
    fn add_in_place(&mut self, other: &Self) {
        self[0] = self[0] + other[0];
        self[1] = self[1] + other[1];
        self[2] = self[2] + other[2];
        self[3] = self[3] + other[3];
    }
    fn sub(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i] - other[i])
    }

    fn mult_mat_col_major(&self, other: &Self) -> Self {
        [
            self[0] * other[0] + self[2] * other[1],
            self[1] * other[0] + self[3] * other[1],
            self[0] * other[2] + self[2] * other[3],
            self[1] * other[2] + self[3] * other[3],
        ]
    }

    fn mult_vec(&self, other: &[Real; 2]) -> [Real; 2] {
        mult_vec(self, other)
    }

    fn transpose(&self) -> Self {
        [self[0], self[2], self[1], self[3]]
    }
    fn squared_norm(&self) -> Real {
        squared_norm(self)
    }
    fn determinant(&self) -> Real {
        determinant(self)
    }

    fn set_zero(&mut self) {
        self.fill(Real::zero());
    }
}

// -------------------------------
// below: from methods

pub fn from_outer_product<T>(a: &[T; 2], b: &[T; 2]) -> [T; 4]
where
    T: num_traits::Float,
{
    [a[0] * b[0], a[1] * b[0], a[0] * b[1], a[1] * b[1]]
}

pub fn from_columns<T>(c0: &[T; 2], c1: &[T; 2]) -> [T; 4]
where
    T: num_traits::Float,
{
    [c0[0], c0[1], c1[0], c1[1]]
}

pub fn from_diagonal<T>(d: &[T; 2]) -> [T; 4]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    [d[0], zero, zero, d[1]]
}

pub fn from_identity<T>() -> [T; 4]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, one]
}

// ----------------------------

pub fn to_columns<T>(a: &[T; 4]) -> ([T; 2], [T; 2])
where
    T: num_traits::Float,
{
    ([a[0], a[1]], [a[2], a[3]])
}

// -----------------------------

pub fn determinant<T>(a: &[T; 4]) -> T
where
    T: num_traits::Float,
{
    a[0] * a[3] - a[1] * a[2]
}

pub fn mult_vec<T>(a: &[T; 4], b: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    [a[0] * b[0] + a[2] * b[1], a[1] * b[0] + a[3] * b[1]]
}

/// Add four 2x2 matrices
pub fn add_four<T>(a: &[T; 4], b: &[T; 4], c: &[T; 4], d: &[T; 4]) -> [T; 4]
where
    T: num_traits::Float,
{
    std::array::from_fn(|i| a[i] + b[i] + c[i] + d[i])
}

pub fn transpose<T>(a: &[T; 4]) -> [T; 4]
where
    T: num_traits::Float,
{
    [a[0], a[2], a[1], a[3]]
}

pub fn squared_norm<T>(a: &[T; 4]) -> T
where
    T: num_traits::Float,
{
    a.iter().map(|&v| v * v).fold(T::zero(), |a, b| a + b)
}

/// Singular Value Decomposition (SVD)
/// input = U * S * V^t
///
/// # Returns
/// (U, S, V)
pub fn svd<Real>(f: &[Real; 4]) -> Option<([Real; 4], [Real; 2], [Real; 4])>
where
    Real: num_traits::Float + std::fmt::Debug,
{
    use crate::vec2::Vec2;
    let zero = Real::zero();
    let ft_f = f.transpose().mult_mat_col_major(f);
    let ft_f = crate::mat2_sym::from_mat2_by_symmetrization(&ft_f);
    let (v, mut lm) = crate::mat2_sym::eigen_decomposition(&ft_f);
    lm.iter_mut().for_each(|x| *x = (*x).max(zero).sqrt());
    let u = f.mult_mat_col_major(&v);
    let (mut u0, mut u1) = to_columns(&u);
    u0.normalize_in_place();
    u1.normalize_in_place();
    let u = from_columns(&u0, &u1);
    Some((u, lm, v))
}

pub fn enforce_rotation_matrix_for_svd<Real>(
    u: &[Real; 4],
    l: &[Real; 2],
    v: &[Real; 4],
) -> ([Real; 4], [Real; 2], [Real; 4])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    if determinant(v) < Real::zero() || determinant(u) < Real::zero() {
        let mut u = *u;
        let mut l = *l;
        let mut v = *v;
        if determinant(&v) < Real::zero() {
            v[2] = -v[2]; // v[0,2] = v[0*3+2]
            v[3] = -v[3]; // v[1,2] = v[1*3+2]
            l[1] = -l[1];
        }
        if determinant(&u) < Real::zero() {
            u[2] = -u[2]; // v[0,2] = v[0*3+2]
            u[3] = -u[3]; // v[1,2] = v[1*3+2]
            l[1] = -l[1];
        }
        (u, l, v)
    } else {
        (*u, *l, *v)
    }
}

#[test]
fn test_svd() {
    use crate::mat2_col_major::Mat2ColMajor;
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0);
    for (_iter, is_rot) in itertools::iproduct!(0..100, 0..2) {
        let m: [f64; 4] = std::array::from_fn(|_| rng.random_range(-1f64..1f64));
        let (u, s, v) = svd(&m).unwrap();
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
            let diff = transpose(&u)
                .mult_mat_col_major(&u)
                .sub(&from_identity())
                .squared_norm();
            assert!(diff < 1.0e-20f64, "{}", diff);
        }
        {
            // test V V^t = I
            let diff = transpose(&v)
                .mult_mat_col_major(&v)
                .sub(&from_identity())
                .squared_norm();
            assert!(diff < 1.0e-20f64, "{}", diff);
        }
        {
            // test A = USV^t
            let s = from_diagonal(&s);
            let m1 = u.mult_mat_col_major(&s).mult_mat_col_major(&transpose(&v));
            let diff = m1.sub(&m).squared_norm();
            assert!(diff < 1.0e-20f64, "{diff} {m:?} {m1:?}");
        }
    }
}
