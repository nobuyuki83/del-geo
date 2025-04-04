//! methods for 3D vector

/// trait for 3D vector
pub trait Vec3<Real>
where
    Self: Sized,
{
    fn normalize(&self) -> Self;
    fn scale(&self, s: Real) -> Self;
    fn scale_in_place(&mut self, s: Real);
    fn norm(&self) -> Real;
    fn squared_norm(&self) -> Real;
    fn dot(&self, other: &Self) -> Real;
    fn sub(&self, other: &Self) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn add_in_place(&mut self, other: &Self);
    fn cross(&self, other: &Self) -> Self;
    fn orthogonalize(&self, v: &Self) -> Self;
    fn transform_homogeneous(&self, v: &[Real; 16]) -> Option<Self>;
    fn xy(&self) -> [Real; 2];
    fn normalize_in_place(&mut self) -> Real;
    fn element_wise_mult(&self, other: &Self) -> Self;
    fn cross_mut(&mut self, v1: &Self, v2: &Self);
}

impl<Real> Vec3<Real> for [Real; 3]
where
    Real: num_traits::Float,
{
    fn normalize(&self) -> Self {
        normalize(self)
    }
    fn scale(&self, s: Real) -> Self {
        scale(self, s)
    }
    fn scale_in_place(&mut self, s: Real) {
        scale_in_place(self, s)
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
    fn add_in_place(&mut self, other: &Self) {
        add_in_place(self, other);
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
    fn normalize_in_place(&mut self) -> Real {
        normalize_in_place(self)
    }
    fn element_wise_mult(&self, other: &Self) -> Self {
        element_wise_mult(self, other)
    }
    fn cross_mut(&mut self, v1: &Self, v2: &Self) {
        cross_mut(self, v1, v2)
    }
}

/// orthogonalize v against u, remove u component from v
/// v - u.scale( dot(u, v) / dot(u, u) )
pub fn orthogonalize<Real>(u: &[Real; 3], v: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let t = u.dot(v) / u.dot(u);
    v.sub(&u.scale(t))
}

/// From axis-angle vector returns a rotation matrix
pub fn to_mat3_from_axisangle_vec<T>(vec: &[T; 3]) -> [T; 9]
where
    T: num_traits::Float,
{
    let one = T::one();
    let sqt = vec.squared_norm();
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
    let aa1 = crate::mat3_col_major::to_vec3_axisangle_from_rot_mat(&rmat);
    let aa0 = del_geo_nalgebra::vec3::from_array(&aa0);
    let aa1 = del_geo_nalgebra::vec3::from_array(&aa1);
    assert!((aa0 - aa1).norm() < 1.0e-5);
}

fn to_mat3_from_dw_normalize<T>(u: &[T; 3]) -> [T; 9]
where
    T: num_traits::Float,
{
    use crate::mat3_col_major::Mat3ColMajor;
    use Vec3;
    let s = T::one() / u.norm();
    let a = crate::mat3_col_major::from_scaled_outer_product(s * s, u, u);
    let b = crate::mat3_col_major::from_identity();
    b.sub(&a).scale(s)
}

// above: To method
// -------------------------------------

/// `vec_n` can be not-normalized vector
pub fn basis_xy_from_basis_z<Real>(vec_n: &[Real; 3]) -> ([Real; 3], [Real; 3])
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let zero = Real::zero();
    let vec_s: [Real; 3] = [zero, one, zero];
    let vec_x = vec_s.cross(vec_n);
    let len = vec_x.norm();
    if len < Real::epsilon() {
        let vec_n = vec_n.normalize();
        let vec_t = [one, zero, zero];
        let vec_x = vec_t.cross(&vec_n); // z????
        let vec_y = vec_n.cross(&vec_x); // x????
        (vec_x, vec_y)
    } else {
        let invlen = one / len;
        let vec_x = vec_x.scale(invlen);
        let vec_y = vec_n.cross(&vec_x);
        (vec_x, vec_y)
    }
}

#[test]
fn test_basis_xy_from_basis_z() {
    let vec_z = [0.3, 0.1, -0.5].normalize();
    let (vec_x, vec_y) = basis_xy_from_basis_z(&vec_z);
    let bases = [vec_x, vec_y, vec_z];
    for i in 0..3 {
        for j in i..3 {
            let dij: f64 = bases[i].dot(&bases[j]);
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

/// axis-aligned basis
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
    std::array::from_fn(|i| a[i] + b[i])
}

pub fn add_in_place<T>(a: &mut [T; 3], b: &[T; 3])
where
    T: num_traits::Float,
{
    *a = a.add(b);
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
    v.squared_norm().sqrt()
}

/// in-place normalize function
pub fn normalize_in_place<T>(v: &mut [T; 3]) -> T
where
    T: num_traits::Float,
{
    let l = v.norm();
    let linv = T::one() / l;
    *v = v.scale(linv);
    l
}

/// return normalized 3D vector
pub fn normalize<T>(v: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    let l = v.norm();
    let linv = T::one() / l;
    v.scale(linv)
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
    std::array::from_fn(|i| a[i] - b[i])
}

pub fn scale<T>(a: &[T; 3], s: T) -> [T; 3]
where
    T: Copy + std::ops::Mul<Output = T>,
{
    std::array::from_fn(|i| a[i] * s)
}

pub fn scale_in_place<T>(a: &mut [T; 3], s: T)
where
    T: num_traits::Float,
{
    *a = a.scale(s);
}

pub fn distance<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    p1.sub(p0).norm()
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
    x.scale(alpha).add(y)
}

pub fn to_quaternion_from_axis_angle_vector<Real>(a: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let half = one / two;
    let one8th = one / (two * two * two);
    let sqlen = a.squared_norm();
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

pub fn mirror_reflection<Real>(v: &[Real; 3], nrm: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let a = nrm.dot(v);
    std::array::from_fn(|i| v[i] - nrm[i] * Real::from(2).unwrap() * a)
}

pub fn element_wise_mult<Real>(a: &[Real; 3], b: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    std::array::from_fn(|i| a[i] * b[i])
}

pub fn add_three<T>(a: &[T; 3], b: &[T; 3], c: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]]
}

pub fn sample_unit_cube<Reng, T>(rng: &mut Reng) -> [T; 3]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<T>,
{
    std::array::from_fn(|_i| rng.random())
}

pub fn wdw_angle_between_two_vecs_using_half_tan<T>(
    u: &[T; 3],
    v: &[T; 3],
) -> ([T; 3], [T; 9], [T; 9])
where
    T: num_traits::Float + std::fmt::Debug,
{
    use crate::mat3_col_major::Mat3ColMajor;
    let uu = u.normalize();
    let uv = v.normalize();
    let c = uu.dot(&uv); // cosine
    let ns = cross(&uu, &uv);
    let b = T::one() + c;
    let w = ns.scale(T::one() / b);
    // dbg!(w, ns, b, bsq_inv);
    let duu_du = to_mat3_from_dw_normalize(u);
    let dvv_dv = to_mat3_from_dw_normalize(v);
    //
    let bsq_inv = T::one() / (b * b);
    let dw_du = {
        let dc_du = duu_du.mult_vec(&uv);
        let mat2 = crate::mat3_col_major::from_vec3_to_skew_mat(&uv);
        let mat0 = mat2.mult_mat_col_major(&duu_du).scale(-b);
        let mat1 = crate::mat3_col_major::from_scaled_outer_product(T::one(), &ns, &dc_du);
        mat0.sub(&mat1).scale(bsq_inv)
    };
    let dw_dv = {
        let dc_dv = dvv_dv.mult_vec(&uu);
        let mat2 = crate::mat3_col_major::from_vec3_to_skew_mat(&uu);
        let mat0 = mat2.mult_mat_col_major(&dvv_dv).scale(b);
        let mat1 = crate::mat3_col_major::from_scaled_outer_product(T::one(), &ns, &dc_dv);
        mat0.sub(&mat1).scale(bsq_inv)
    };
    (w, dw_du, dw_dv)
}

#[test]
fn test_wdw_angle_between_two_vecs_using_half_tan() {
    use Vec3;
    let u0 = [0.3, -0.5, 0.1];
    let v0 = [-0.4, 0.2, -0.7];
    let (w0, dw0du, dw0dv) = wdw_angle_between_two_vecs_using_half_tan(&u0, &v0);
    dbg!(w0);
    dbg!(dw0du);
    dbg!(dw0dv);
    let eps = 1.0e-5f64;
    for i_dim in 0..3 {
        {
            let u1 = {
                let mut u1 = u0;
                u1[i_dim] += eps;
                u1
            };
            let (w1, _dw1du, _dw1dv) = wdw_angle_between_two_vecs_using_half_tan(&u1, &v0);
            let diff_num = w1.sub(&w0).scale(1.0 / eps);
            let diff_ana = crate::mat3_col_major::to_vec3_column(&dw0du, i_dim);
            let diffnrm = diff_ana.sub(&diff_num).norm();
            assert!(
                diffnrm < diff_ana.norm() * 1.0e-4,
                "{} {}",
                diffnrm,
                diff_ana.norm()
            );
        }
        //
        {
            let v1 = {
                let mut v1 = v0;
                v1[i_dim] += eps;
                v1
            };
            let (w1, _dw1du, _dw1dv) = wdw_angle_between_two_vecs_using_half_tan(&u0, &v1);
            let diff_num = w1.sub(&w0).scale(1.0 / eps);
            let diff_ana = crate::mat3_col_major::to_vec3_column(&dw0dv, i_dim);
            let diffnrm = diff_ana.sub(&diff_num).norm();
            assert!(
                diffnrm < diff_ana.norm() * 1.0e-4,
                "{} {}",
                diffnrm,
                diff_ana.norm()
            );
        }
    }
}

// ------------------------------------------
#[derive(Debug, Clone, Copy)]
pub struct XYZ<'a, Real> {
    pub p: &'a [Real; 3],
}

impl<Real> XYZ<'_, Real>
where
    Real: num_traits::Float,
{
    pub fn aabb(&self) -> [Real; 6] {
        [
            self.p[0], self.p[1], self.p[2], self.p[0], self.p[1], self.p[2],
        ]
    }
}
