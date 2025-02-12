//! methods for 2D vector

/// trait for 2D vector
pub trait Vec2<Real>
where
    Self: Sized,
{
    fn sub(&self, other: &Self) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn add_in_place(&mut self, other: &Self);
    fn transform_homogeneous(&self, v: &[Real; 9]) -> Option<[Real; 2]>;
    fn dot(&self, other: &Self) -> Real;
    fn scale(&self, s: Real) -> Self;
    fn scale_in_place(&mut self, s: Real);
    fn orthogonalize(&self, v: &Self) -> Self;
    fn norm(&self) -> Real;
    fn squared_norm(&self) -> Real;
    fn aabb(&self) -> [Real; 4];
    fn normalize(&self) -> Self;
    fn normalize_in_place(&mut self);
    fn cross(&self, other: &Self) -> Real;
}

impl<Real> Vec2<Real> for [Real; 2]
where
    Real: num_traits::Float,
{
    fn sub(&self, other: &Self) -> Self {
        sub(self, other)
    }
    fn add(&self, other: &Self) -> Self {
        add(self, other)
    }
    fn add_in_place(&mut self, other: &Self) {
        add_in_place(self, other);
    }
    fn transform_homogeneous(&self, v: &[Real; 9]) -> Option<Self> {
        crate::mat3_col_major::transform_homogeneous(v, self)
    }
    fn dot(&self, other: &Self) -> Real {
        dot(self, other)
    }
    fn scale(&self, s: Real) -> Self {
        scale(self, s)
    }
    fn scale_in_place(&mut self, s: Real) {
        scale_in_place(self, s);
    }
    fn orthogonalize(&self, v: &Self) -> Self {
        orthogonalize(self, v)
    }
    fn norm(&self) -> Real {
        length(self)
    }
    fn squared_norm(&self) -> Real {
        squared_length(self)
    }
    fn aabb(&self) -> [Real; 4] {
        aabb(self)
    }
    fn normalize(&self) -> Self {
        normalize(self)
    }
    fn normalize_in_place(&mut self) {
        normalize_in_place(self);
    }
    /// the z coordiante of the cross product of two vectors in the xy-plane
    fn cross(&self, other: &Self) -> Real {
        cross(self, other)
    }
}

pub fn basis<T>(i_dim: usize, eps: T) -> [T; 2]
where
    T: num_traits::Float,
{
    let mut b = [T::zero(); 2];
    b[i_dim] = eps;
    b
}

pub fn length<Real>(p: &[Real; 2]) -> Real
where
    Real: num_traits::Float,
{
    (p[0] * p[0] + p[1] * p[1]).sqrt()
}

pub fn squared_length<Real>(p: &[Real; 2]) -> Real
where
    Real: num_traits::Float,
{
    p[0] * p[0] + p[1] * p[1]
}

pub fn sub<T>(a: &[T; 2], b: &[T; 2]) -> [T; 2]
where
    T: std::ops::Sub<Output = T> + Copy,
{
    std::array::from_fn(|i| a[i] - b[i])
}

pub fn add<T>(a: &[T; 2], b: &[T; 2]) -> [T; 2]
where
    T: std::ops::Add<Output = T> + Copy,
{
    std::array::from_fn(|i| a[i] + b[i])
}

pub fn add_in_place<T>(a: &mut [T; 2], b: &[T; 2])
where
    T: std::ops::Add<Output = T> + Copy,
{
    a.iter_mut().zip(b.iter()).for_each(|(a, &b)| *a = *a + b);
}

pub fn rotate90<T>(v: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float + Copy,
{
    [-v[1], v[0]]
}

pub fn add_three<T>(a: &[T; 2], b: &[T; 2], c: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    [a[0] + b[0] + c[0], a[1] + b[1] + c[1]]
}

pub fn scale<T>(a: &[T; 2], s: T) -> [T; 2]
where
    T: num_traits::Float,
{
    std::array::from_fn(|i| a[i] * s)
}

pub fn scale_in_place<T>(a: &mut [T; 2], s: T)
where
    T: num_traits::Float,
{
    a.iter_mut().for_each(|v| *v = *v * s);
}

pub fn dot<T>(a: &[T; 2], b: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    a[0] * b[0] + a[1] * b[1]
}

pub fn area_quadrilateral<T>(a: &[T; 2], b: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    a[0] * b[1] - a[1] * b[0]
}

pub fn angle_between_two_vecs<T>(a: &[T; 2], b: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    let dot = a.dot(b);
    let area = area_quadrilateral(a, b);
    area.atan2(dot)
}

#[test]
fn test_angle_between_two_vecs() {
    let a = [3f64.sqrt(), 1.0];
    let b = [-1.0, 1.0];
    let theta0 = angle_between_two_vecs(&a, &b);
    let theta1 = 7f64 / 12f64 * std::f64::consts::PI;
    assert!((theta0 - theta1).abs() < 1.0e-10);
}

pub fn wdw_angle_between_two_vecs<T>(u: &[T; 2], v: &[T; 2]) -> (T, [[T; 2]; 2])
where
    T: num_traits::Float,
{
    let a = dot(u, v);
    let b = area_quadrilateral(u, v);
    let w = b.atan2(a);
    let tmp0 = T::one() / (a * a + b * b);
    let dw_da = -b * tmp0;
    let dw_db = a * tmp0;
    let dw_du = [dw_da * v[0] + dw_db * v[1], dw_da * v[1] - dw_db * v[0]];
    let dw_dv = [dw_da * u[0] - dw_db * u[1], dw_da * u[1] + dw_db * u[0]];
    (w, [dw_du, dw_dv])
}

#[test]
fn test_wdw_angle_between_two_vecs() {
    let a0 = [[3f64.sqrt(), 1.0], [-1.0, 1.0]];
    let (t0, dt0) = wdw_angle_between_two_vecs(&a0[0], &a0[1]);
    let eps = 1.0e-5;
    for (ino, idim) in itertools::iproduct!(0..2, 0..2) {
        let a1 = {
            let mut a1 = a0;
            a1[ino][idim] += eps;
            a1
        };
        let (t1, _dt1) = wdw_angle_between_two_vecs(&a1[0], &a1[1]);
        let v0 = (t1 - t0) / eps;
        let v1 = dt0[ino][idim];
        assert!((v0 - v1).abs() < 1.0e-5);
    }
}

/// convert homogeneous coordinate to the Cartesian coordiante
pub fn from_homogeneous<Real>(v: &[Real; 3]) -> Option<[Real; 2]>
where
    Real: num_traits::Float,
{
    if v[2].is_zero() {
        return None;
    }
    Some([v[0] / v[2], v[0] / v[2]])
}

/// rotate around the origin
/// # Argument
/// - theta: radian
pub fn rotate<Real>(p: &[Real; 2], theta: Real) -> [Real; 2]
where
    Real: num_traits::Float,
{
    let c = theta.cos();
    let s = theta.sin();
    [c * p[0] - s * p[1], s * p[0] + c * p[1]]
}

pub fn normalize<Real>(p: &[Real; 2]) -> [Real; 2]
where
    Real: num_traits::Float,
{
    let invl = Real::one() / (p[0] * p[0] + p[1] * p[1]).sqrt();
    p.scale(invl)
}

pub fn normalize_in_place<Real>(p: &mut [Real; 2])
where
    Real: num_traits::Float,
{
    let invl = Real::one() / (p[0] * p[0] + p[1] * p[1]).sqrt();
    p.scale_in_place(invl);
}

pub fn orthogonalize<T>(u: &[T; 2], v: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    let t = u.dot(v) / u.dot(u);
    v.sub(&u.scale(t))
}

pub fn axpy<Real>(alpha: Real, x: &[Real; 2], y: &[Real; 2]) -> [Real; 2]
where
    Real: num_traits::Float,
{
    x.scale(alpha).add(y)
}

pub fn aabb<T>(p: &[T; 2]) -> [T; 4]
where
    T: num_traits::Float,
{
    [p[0], p[1], p[0], p[1]]
}

/// the z coordinate of the cross product of two vectors in the xy-plane
pub fn cross<T>(a: &[T; 2], b: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    a[0] * b[1] - b[0] * a[1]
}

// -------------------------------
// below: about the Vec2 class
#[derive(Debug, Clone, Copy)]
pub struct XY<Real> {
    pub p: [Real; 2],
}

impl<Real> XY<Real>
where
    Real: num_traits::Float,
{
    pub fn new(p: [Real; 2]) -> Self {
        Self { p }
    }
}

impl<Real> std::ops::Add<XY<Real>> for XY<Real>
where
    Real: num_traits::Float,
{
    type Output = XY<Real>;
    fn add(self, rhs: Self) -> XY<Real> {
        XY {
            p: std::array::from_fn(|i| self.p[i] + rhs.p[i]),
        }
    }
}

impl<Real> std::ops::Sub<XY<Real>> for XY<Real>
where
    Real: num_traits::Float,
{
    type Output = XY<Real>;
    fn sub(self, rhs: Self) -> XY<Real> {
        XY {
            p: std::array::from_fn(|i| self.p[i] - rhs.p[i]),
        }
    }
}
