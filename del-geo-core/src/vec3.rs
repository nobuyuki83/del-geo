//! methods for 3D vector
//!
use std::ops::MulAssign;
pub trait Vec3<Real>
where
    Self: Sized,
{
    fn normalized(&self) -> Self;
    fn scaled(&self, s: Real) -> Self;
    fn norm(&self) -> Real;
    fn squared_norm(&self) -> Real;
    fn dot(&self, other: &Self) -> Real;
    fn sub(&self, other: &Self) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn cross(&self, other: &Self) -> Self;
    fn orthogonalize(&self, v: &Self) -> Self;
    fn transform_homogeneous(&self, v: &[Real; 16]) -> Option<Self>;
    fn xy(&self) -> [Real; 2];
    fn normalize(&mut self) -> Real;
}

impl<Real> Vec3<Real> for [Real; 3]
where
    Real: num_traits::Float + std::ops::MulAssign,
{
    fn normalized(&self) -> Self {
        normalized(self)
    }
    fn scaled(&self, s: Real) -> Self {
        scaled(self, s)
    }
    fn squared_norm(&self) -> Real {
        squared_norm(self)
    }
    fn sub(&self, other: &Self) -> Self {
        sub(self, other)
    }
    fn add(&self, other: &Self) -> Self {
        add(self, other)
    }
    fn dot(&self, other: &Self) -> Real {
        dot(self, other)
    }
    fn norm(&self) -> Real {
        norm(self)
    }
    fn cross(&self, other: &Self) -> Self {
        cross(self, other)
    }
    fn orthogonalize(&self, v: &Self) -> Self {
        orthogonalize(self, v)
    }
    fn transform_homogeneous(&self, v: &[Real; 16]) -> Option<[Real; 3]> {
        crate::mat4_col_major::transform_homogeneous(v, self)
    }
    fn xy(&self) -> [Real; 2] {
        [self[0], self[1]]
    }
    fn normalize(&mut self) -> Real {
        normalize(self)
    }
}

pub fn orthogonalize<Real>(u: &[Real; 3], v: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float + MulAssign,
{
    let t = u.dot(v) / u.dot(u);
    [v[0] - t * u[0], v[1] - t * u[1], v[2] - t * u[2]]
}

pub fn to_mat3_from_axisangle_vec<T>(vec: &[T; 3]) -> [T; 9]
where
    T: num_traits::Float,
{
    let one = T::one();
    let sqt = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    if sqt <= T::epsilon() {
        // infinitesimal rotation approximation
        return [
            one, vec[2], -vec[1], -vec[2], one, vec[0], vec[1], -vec[0], one,
        ];
    }
    let t = sqt.sqrt();
    let invt = one / t;
    let n = [vec[0] * invt, vec[1] * invt, vec[2] * invt];
    let c0 = t.cos();
    let s0 = t.sin();
    [
        c0 + (one - c0) * n[0] * n[0],
        n[2] * s0 + (one - c0) * n[1] * n[0],
        -n[1] * s0 + (one - c0) * n[2] * n[0],
        -n[2] * s0 + (one - c0) * n[0] * n[1],
        c0 + (one - c0) * n[1] * n[1],
        n[0] * s0 + (one - c0) * n[2] * n[1],
        n[1] * s0 + (one - c0) * n[0] * n[2],
        -n[0] * s0 + (one - c0) * n[1] * n[2],
        c0 + (one - c0) * n[2] * n[2],
    ]
}

#[test]
fn test_to_mat3_from_axisangle_vec() {
    let aa0 = [1.0, 0.3, -0.5f64];
    let rmat = to_mat3_from_axisangle_vec(&aa0);
    let aa1 = crate::mat3_col_major::to_vec3_axisangle_from_rot_mat(rmat);
    let aa0 = del_geo_nalgebra::vec3::from_array(&aa0);
    let aa1 = del_geo_nalgebra::vec3::from_array(&aa1);
    assert!((aa0 - aa1).norm() < 1.0e-5);
}

pub fn basis_xy_from_basis_z<Real>(vec_n: &[Real; 3]) -> ([Real; 3], [Real; 3])
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let zero = Real::zero();
    let vec_s: [Real; 3] = [zero, one, zero];
    let vec_x = cross(&vec_s, vec_n);
    let len = norm(&vec_x);
    if len < Real::epsilon() {
        let vec_t = [one, zero, zero];
        let vec_x = cross(&vec_t, vec_n); // z????
        let vec_y = cross(vec_n, &vec_x); // x????
        (vec_x, vec_y)
    } else {
        let invlen = one / len;
        let vec_x = [vec_x[0] * invlen, vec_x[1] * invlen, vec_x[2] * invlen];
        let vec_y = cross(vec_n, &vec_x);
        (vec_x, vec_y)
    }
}

#[test]
fn test_basis_xy_from_basis_z() {
    let vec_z = normalized(&[0.3, 0.1, -0.5]);
    let (vec_x, vec_y) = basis_xy_from_basis_z(&vec_z);
    let bases = [vec_x, vec_y, vec_z];
    for i in 0..3 {
        for j in i..3 {
            let dij: f64 = dot(&bases[i], &bases[j]);
            if i == j {
                assert!((dij - 1.0).abs() < 1.0e-5);
            } else {
                assert!(dij.abs() < 1.0e-5);
            }
        }
    }
    let stp = scalar_triple_product(&vec_x, &vec_y, &vec_z);
    assert!((stp - 1.0).abs() < 1.0e-5);
}

pub fn basis<Real>(i_dim: usize, eps: Real) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    match i_dim {
        0 => [eps, zero, zero],
        1 => [zero, eps, zero],
        2 => [zero, zero, eps],
        _ => panic!(),
    }
}

pub fn add<T>(a: &[T; 3], b: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn squared_norm<T>(p: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    assert_eq!(p.len(), 3);
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
}

pub fn norm<T>(v: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// in-place normalize function
pub fn normalize<T>(v: &mut [T; 3]) -> T
where
    T: num_traits::Float + std::ops::MulAssign,
{
    let l = norm(v);
    let linv = T::one() / l;
    v[0] *= linv;
    v[1] *= linv;
    v[2] *= linv;
    l
}

/// in-place normalize function
pub fn normalized<T>(v: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    let l = norm(v);
    let linv = T::one() / l;
    std::array::from_fn(|i| v[i] * linv)
}

pub fn cross_mut<T>(vo: &mut [T; 3], v1: &[T; 3], v2: &[T; 3])
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + Copy,
{
    vo[0] = v1[1] * v2[2] - v2[1] * v1[2];
    vo[1] = v1[2] * v2[0] - v2[2] * v1[0];
    vo[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

pub fn cross<T>(v1: &[T; 3], v2: &[T; 3]) -> [T; 3]
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + Copy,
{
    [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1],
    ]
}

pub fn dot<T>(a: &[T; 3], b: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// subtracting two 3d vectors
/// * `a` - 3d vector
/// * `b` - 3d vector
///   return a-b
pub fn sub<T>(a: &[T; 3], b: &[T; 3]) -> [T; 3]
where
    T: std::ops::Sub<Output = T> + Copy,
{
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scaled<T>(a: &[T; 3], s: T) -> [T; 3]
where
    T: Copy + std::ops::Mul<Output = T>,
{
    [s * a[0], s * a[1], s * a[2]]
}

pub fn scale<T>(a: &mut [T; 3], s: T)
where
    T: MulAssign + Copy,
{
    for x in a.iter_mut() {
        *x *= s;
    }
}

pub fn distance<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let v0 = p1[0] - p0[0];
    let v1 = p1[1] - p0[1];
    let v2 = p1[2] - p0[2];
    (v0 * v0 + v1 * v1 + v2 * v2).sqrt()
}

pub fn scalar_triple_product<T>(a: &[T; 3], b: &[T; 3], c: &[T; 3]) -> T
where
    T: std::ops::Mul<Output = T> + std::ops::Sub<Output = T> + std::ops::Add<Output = T> + Copy,
{
    let v0: T = a[0] * (b[1] * c[2] - b[2] * c[1]);
    let v1: T = a[1] * (b[2] * c[0] - b[0] * c[2]);
    let v2: T = a[2] * (b[0] * c[1] - b[1] * c[0]);
    v0 + v1 + v2
}

pub fn axpy<Real>(alpha: Real, x: &[Real; 3], y: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    [
        alpha * x[0] + y[0],
        alpha * x[1] + y[1],
        alpha * x[2] + y[2],
    ]
}

pub fn to_quaternion_from_axis_angle_vector<Real>(a: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let half = one / two;
    let one8th = one / (two * two * two);
    let sqlen = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if sqlen <= Real::epsilon() {
        return [half * a[0], half * a[1], half * a[2], one - one8th * sqlen];
    }
    let lena = sqlen.sqrt();
    [
        (lena * half).sin() * a[0] / lena,
        (lena * half).sin() * a[1] / lena,
        (lena * half).sin() * a[2] / lena,
        (lena * half).cos(),
    ]
}

pub fn mirror_reflection(v: &[f32; 3], nrm: &[f32; 3]) -> [f32; 3] {
    let a = dot(nrm, v);
    [
        v[0] - nrm[0] * 2. * a,
        v[1] - nrm[1] * 2. * a,
        v[2] - nrm[2] * 2. * a,
    ]
}

pub fn element_wise_mult(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

// ------------------------------------------
pub struct XYZ<'a, Real> {
    pub p: &'a [Real; 3],
}

#[allow(clippy::needless_lifetimes)]
impl<'a, Real> XYZ<'a, Real>
where
    Real: num_traits::Float,
{
    pub fn aabb(&self) -> [Real; 6] {
        [
            self.p[0], self.p[1], self.p[2], self.p[0], self.p[1], self.p[2],
        ]
    }
}
