//! methods for 4x4 matrix

use std::ops::AddAssign;

use num_traits::AsPrimitive;

pub fn transform_homogeneous<Real>(transform: &[Real; 16], x: &[Real; 3]) -> Option<[Real; 3]>
    where
        Real: num_traits::Float,
{
    let y3 = transform[3] * x[0] + transform[7] * x[1] + transform[11] * x[2] + transform[15];
    if y3.is_zero() {
        return None;
    }
    //
    let y0 = transform[0] * x[0] + transform[4] * x[1] + transform[8] * x[2] + transform[12];
    let y1 = transform[1] * x[0] + transform[5] * x[1] + transform[9] * x[2] + transform[13];
    let y2 = transform[2] * x[0] + transform[6] * x[1] + transform[10] * x[2] + transform[14];
    Some([y0 / y3, y1 / y3, y2 / y3])
}

pub fn transform_vector<Real>(transform: &[Real; 16], x: &[Real; 3]) -> [Real; 3]
    where
        Real: num_traits::Float,
{
    let y0 = transform[0] * x[0] + transform[4] * x[1] + transform[8] * x[2];
    let y1 = transform[1] * x[0] + transform[5] * x[1] + transform[9] * x[2];
    let y2 = transform[2] * x[0] + transform[6] * x[1] + transform[10] * x[2];
    [y0, y1, y2]
}

pub fn identity<Real>() -> [Real; 16]
    where
        Real: num_traits::Zero + num_traits::One + Copy,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        one, zero, zero, zero, zero, one, zero, zero, zero, zero, one, zero, zero, zero, zero, one,
    ]
}

pub fn diagonal<Real>(m11: Real, m22: Real, m33: Real, m44: Real) -> [Real; 16]
    where
        Real: num_traits::Zero + Copy,
{
    let zero = Real::zero();
    [
        m11, zero, zero, zero, zero, m22, zero, zero, zero, zero, m33, zero, zero, zero, zero, m44,
    ]
}

pub fn try_inverse<Real>(b: &[Real; 16]) -> Option<[Real; 16]>
    where
        Real: num_traits::Float + std::ops::MulAssign + std::ops::SubAssign,
{
    crate::matn::try_inverse::<Real, 4, 16>(b)
}

/// perspective transformation matrix (column major) compatible with blender
/// * asp - aspect ratio (width / height)
/// * lens - the forcus distance (unit: mm) where the sensor size for longet edge is 18*2 mm.
/// * near - distancd to the near clipping planne (>0)
/// * far - distance  to the far clipping plane (>0)
pub fn camera_perspective_blender<Real>(asp: Real, lens: Real, near: Real, far: Real) -> [Real; 16]
    where
        Real: num_traits::Float + 'static,
        f64: AsPrimitive<Real>,
{
    let zero = Real::zero();
    let one = Real::one();
    let (a, b) = if asp < one {
        (one / asp, one)
    } else {
        (one, asp)
    };
    let half_theta = 18f64.as_().atan2(lens);
    let focus = half_theta.cos() / half_theta.sin();
    let tmp0 = (near + far) / (near - far);
    let tmp1 = (one + one) * far * near / (near - far);
    [
        -a * focus,
        zero,
        zero,
        zero,
        zero,
        -b * focus,
        zero,
        zero,
        zero,
        zero,
        tmp0,
        one,
        zero,
        zero,
        tmp1,
        zero,
    ]
}

pub fn camera_external_blender<Real>(
    cam_location: &[Real; 3],
    cam_rot_x_deg: Real,
    cam_rot_y_deg: Real,
    cam_rot_z_deg: Real,
) -> [Real; 16]
    where
        Real: num_traits::Float + num_traits::FloatConst + 'static + AddAssign,
        f64: AsPrimitive<Real>,
{
    let deg2rad: Real = Real::PI() / 180.0.as_();
    let transl =
        crate::mat4_col_major::translate(&[-cam_location[0], -cam_location[1], -cam_location[2]]);
    let rot_x = crate::mat4_col_major::rot_x(-cam_rot_x_deg * deg2rad);
    let rot_y = crate::mat4_col_major::rot_y(-cam_rot_y_deg * deg2rad);
    let rot_z = crate::mat4_col_major::rot_z(-cam_rot_z_deg * deg2rad);
    let rot_yz = crate::mat4_col_major::multmat(&rot_y, &rot_z);
    let rot_zyx = crate::mat4_col_major::multmat(&rot_x, &rot_yz);
    crate::mat4_col_major::multmat(&rot_zyx, &transl)
}

pub fn translate<Real>(v: &[Real; 3]) -> [Real; 16]
    where
        Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        one, zero, zero, zero, zero, one, zero, zero, zero, zero, one, zero, v[0], v[1], v[2], one,
    ]
}

pub fn multmat<Real>(a: &[Real; 16], b: &[Real; 16]) -> [Real; 16]
    where
        Real: num_traits::Float + std::ops::AddAssign,
{
    let mut o = [Real::zero(); 16];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                o[i + j * 4] += a[i + k * 4] * b[k + j * 4];
            }
        }
    }
    o
}

#[test]
fn test_inverse_multmat() {
    let a: [f64; 16] = [
        1., 3., 4., 8.,
        3., 5., 5., 2.,
        5., 7., 8., 9.,
        8., 4., 5., 0.];
    let ainv = try_inverse(&a).unwrap();
    let ainv_a = multmat(&ainv, &a);
    for i in 0..4 {
        for j in 0..4 {
            if i == j {
                assert!((1.0 - ainv_a[i * 4 + j]).abs() < 1.0e-5f64);
            } else {
                assert!(ainv_a[i * 4 + j].abs() < 1.0e-5f64);
            }
        }
    }
}

pub fn rot_x<Real>(theta: Real) -> [Real; 16]
    where
        Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let c = theta.cos();
    let s = theta.sin();
    [
        one, zero, zero, zero, zero, c, s, zero, zero, -s, c, zero, zero, zero, zero, one,
    ]
}

pub fn rot_y<Real>(theta: Real) -> [Real; 16]
    where
        Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let c = theta.cos();
    let s = theta.sin();
    [
        c, zero, -s, zero, zero, one, zero, zero, s, zero, c, zero, zero, zero, zero, one,
    ]
}

pub fn rot_z<Real>(theta: Real) -> [Real; 16]
    where
        Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let c = theta.cos();
    let s = theta.sin();
    [
        c, s, zero, zero, -s, c, zero, zero, zero, zero, one, zero, zero, zero, zero, one,
    ]
}

