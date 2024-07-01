//! methods for 3D vector

pub fn squared_norm_<T>(p: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    assert_eq!(p.len(), 3);
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
}

pub fn norm_<T>(v: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn normalize_<T>(v: &mut [T; 3]) -> T
where
    T: num_traits::Float + std::ops::MulAssign,
{
    let l = norm_(v);
    let linv = T::one() / l;
    v[0] *= linv;
    v[1] *= linv;
    v[2] *= linv;
    l
}

pub fn cross_mut_<T>(vo: &mut [T; 3], v1: &[T; 3], v2: &[T; 3])
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + Copy,
{
    vo[0] = v1[1] * v2[2] - v2[1] * v1[2];
    vo[1] = v1[2] * v2[0] - v2[2] * v1[0];
    vo[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

pub fn cross_<T>(v1: &[T; 3], v2: &[T; 3]) -> [T; 3]
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + Copy,
{
    [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1],
    ]
}

pub fn dot_<T>(a: &[T; 3], b: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// substracting two 3d vectors
/// * `a` - 3d vector
/// * `b` - 3d vector
/// return a-b
pub fn sub_<T>(a: &[T; 3], b: &[T; 3]) -> [T; 3]
where
    T: std::ops::Sub<Output = T> + Copy,
{
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale_<T>(a: &[T; 3], s: T) -> [T; 3]
where
    T: Copy + std::ops::Mul<Output = T>,
{
    [s * a[0], s * a[1], s * a[2]]
}

pub fn distance_<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let v0 = p1[0] - p0[0];
    let v1 = p1[1] - p0[1];
    let v2 = p1[2] - p0[2];
    (v0 * v0 + v1 * v1 + v2 * v2).sqrt()
}

pub fn scalar_triple_product_<T>(a: &[T; 3], b: &[T; 3], c: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + std::ops::Add<Output = T> + Copy,
{
    let v0: T = a[0] * (b[1] * c[2] - b[2] * c[1]);
    let v1: T = a[1] * (b[2] * c[0] - b[0] * c[2]);
    let v2: T = a[2] * (b[0] * c[1] - b[1] * c[0]);
    v0 + v1 + v2
}

pub fn axpy_<Real>(alpha: Real, x: &[Real; 3], y: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    [
        alpha * x[0] + y[0],
        alpha * x[1] + y[1],
        alpha * x[2] + y[2],
    ]
}
