// quaternion defined as [i,j,k,w] where w is the real part

use crate::vec3::Vec3;
/// trait for quaternion
pub trait Quaternion<Real>
where
    Self: Sized,
{
    fn to_mat3_col_major(&self) -> [Real; 9];
    fn to_mat4_col_major(&self) -> [Real; 16];
    fn normalized(&self) -> Self;
    fn inverse(&self) -> Self;
    fn mult_quaternion(&self, q: &Self) -> Self;
    fn identity() -> Self;
}
impl<Real> Quaternion<Real> for [Real; 4]
where
    Real: num_traits::Float,
{
    fn to_mat3_col_major(&self) -> [Real; 9] {
        to_mat3_col_major(self)
    }
    fn to_mat4_col_major(&self) -> [Real; 16] {
        to_mat4_col_major(self)
    }
    fn normalized(&self) -> Self {
        normalized(self)
    }
    fn inverse(&self) -> Self {
        inverse(*self)
    }
    fn mult_quaternion(&self, q: &Self) -> Self {
        mult_quaternion(self, q)
    }
    fn identity() -> Self {
        identity()
    }
}
pub fn to_mat3_col_major<Real>(q: &[Real; 4]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let x2 = q[0] * q[0] * two;
    let y2 = q[1] * q[1] * two;
    let z2 = q[2] * q[2] * two;
    let xy = q[0] * q[1] * two;
    let yz = q[1] * q[2] * two;
    let zx = q[2] * q[0] * two;
    let xw = q[0] * q[3] * two;
    let yw = q[1] * q[3] * two;
    let zw = q[2] * q[3] * two;
    [
        one - y2 - z2,
        xy + zw,
        zx - yw,
        xy - zw,
        one - z2 - x2,
        yz + xw,
        zx + yw,
        yz - xw,
        one - x2 - y2,
    ]
}

pub fn to_mat4_col_major<Real>(q: &[Real; 4]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let a = to_mat3_col_major(q);
    [
        a[0], a[1], a[2], zero, a[3], a[4], a[5], zero, a[6], a[7], a[8], zero, zero, zero, zero,
        one,
    ]
}

pub fn normalized<Real>(q: &[Real; 4]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let invlen = Real::one() / len;
    [q[0] * invlen, q[1] * invlen, q[2] * invlen, q[3] * invlen]
}

pub fn inverse<Real>(q: [Real; 4]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let sqlen = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    let sqleninv = Real::one() / sqlen;
    [
        -q[0] * sqleninv,
        -q[1] * sqleninv,
        -q[2] * sqleninv,
        q[3] * sqleninv,
    ]
}

#[test]
fn hoge() {
    use crate::mat3_col_major::Mat3ColMajor;
    // use crate::vec3::Vec3;
    let quat: [f64; 4] = [1., 2., 3., 1.];
    let quat = quat.normalized();
    let r_mat = quat.to_mat3_col_major();
    let r_mat_transp = r_mat.transpose();
    {
        let identity0 = r_mat.mult_mat_col_major(&r_mat_transp);
        // dbg!(&identity0);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert!((identity0[i + 3 * j] - 1.0).abs() < 1.0e-5);
                } else {
                    assert!(
                        identity0[i * 3 + j].abs() < 1.0e-5,
                        "{}",
                        identity0[i * 3 * j]
                    );
                }
            }
        }
    }
    // for 3D Gaussian splatting
    let s_mat = crate::mat3_col_major::from_diagonal(&[1.0f64, 2.0, 0.3]);
    let rs = r_mat.mult_mat_col_major(&s_mat);
    let in_vec = [0.2, 2.3, 0.1];
    let out_vec = rs.transpose().mult_vec(&in_vec);
    let _l = out_vec.squared_norm();
    // dbg!(&out_vec, l);
}

pub fn mult_quaternion<Real>(p: &[Real; 4], q: &[Real; 4]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    [
        p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1],
        p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0],
        p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3],
        p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2],
    ]
}

pub fn from_axisangle<Real>(a: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let half = Real::from(0.5).unwrap();
    let sqlen = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    if sqlen < Real::epsilon() {
        return [
            half * a[0],
            half * a[1],
            half * a[2],
            -sqlen * Real::from(0.125).unwrap() + Real::one(),
        ];
    }
    let lena = sqlen.sqrt();
    [
        (lena * half).sin() * a[0] / lena,
        (lena * half).sin() * a[1] / lena,
        (lena * half).sin() * a[2] / lena,
        (lena * half).cos(),
    ]
}

pub fn identity<Real>() -> [Real; 4]
where
    Real: num_traits::Float,
{
    [Real::zero(), Real::zero(), Real::zero(), Real::one()]
}

/// return rotation around axis with radian
pub fn around_axis<Real>(a: &[Real; 3], rad: Real) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let v = a.normalize();
    let half = rad / Real::from(2).unwrap();
    let sin = half.sin();
    [v[0] * sin, v[1] * sin, v[2] * sin, half.cos()]
}
