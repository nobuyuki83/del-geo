//! methods for 4x4 matrix

use num_traits::AsPrimitive;
use std::ops::AddAssign;

pub fn from_identity<Real>() -> [Real; 16]
where
    Real: num_traits::Zero + num_traits::One + Copy,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        one, zero, zero, zero, zero, one, zero, zero, zero, zero, one, zero, zero, zero, zero, one,
    ]
}

pub fn from_diagonal<Real>(m11: Real, m22: Real, m33: Real, m44: Real) -> [Real; 16]
where
    Real: num_traits::Zero + Copy,
{
    let zero = Real::zero();
    [
        m11, zero, zero, zero, zero, m22, zero, zero, zero, zero, m33, zero, zero, zero, zero, m44,
    ]
}

pub fn from_scale_uniform<Real>(s: Real) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        s, zero, zero, zero, zero, s, zero, zero, zero, zero, s, zero, zero, zero, zero, one,
    ]
}

pub fn from_translate<Real>(v: &[Real; 3]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        one, zero, zero, zero, zero, one, zero, zero, zero, zero, one, zero, v[0], v[1], v[2], one,
    ]
}

pub fn from_rot_x<Real>(theta: Real) -> [Real; 16]
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

pub fn from_rot_y<Real>(theta: Real) -> [Real; 16]
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

pub fn from_rot_z<Real>(theta: Real) -> [Real; 16]
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

/// transformation converting normalized device coordinate (NDC) `[-1,+1]^2` to pixel coordinate
/// depth (-1, +1) is transformed to (0, +1)
/// for example:
///     [-1,-1,-1] becomes (0, H, 0)
///     [+1,+1,+1] becomes (W, 0, 1)
/// * Return
///     * `image_shape` - (width, height)
pub fn from_transform_ndc2pix(img_shape: (usize, usize)) -> [f32; 16] {
    [
        0.5 * (img_shape.0 as f32),
        0.,
        0.,
        0.,
        //
        0.,
        -0.5 * (img_shape.1 as f32),
        0.,
        0.,
        //
        0.,
        0.,
        0.5,
        0.,
        //
        0.5 * (img_shape.0 as f32),
        0.5 * (img_shape.1 as f32),
        0.5,
        1.,
    ]
}

pub fn from_aabb3_fit_into_ndc_preserving_xyasp(aabb: &[f32; 6], asp: f32) -> [f32; 16] {
    let cntr = crate::aabb3::center(aabb);
    let (scale_xy, scale_z) = {
        let size = crate::aabb3::size(aabb);
        let lenx = size[0] / asp;
        let leny = size[1];
        dbg!(size, lenx, leny);
        if lenx > leny {
            (lenx / 2f32, size[2] / 2f32)
        } else {
            (leny / 2f32, size[2] / 2f32)
        }
    };
    [
        scale_xy * asp,
        0.,
        0.,
        0.,
        0.,
        scale_xy,
        0.,
        0.,
        0.,
        0.,
        scale_z,
        0.,
        cntr[0],
        cntr[1],
        cntr[2],
        1.,
    ]
}

// above: from method (making 4x4 matrix)
// ----------------------------------------

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

pub fn jacobian_transform<Real>(t: &[Real; 16], p: &[Real; 3]) -> [Real; 9]
where
    Real: num_traits::Float + Copy,
{
    let a = [t[0], t[1], t[2], t[4], t[5], t[6], t[8], t[9], t[10]];
    let b = [t[12], t[13], t[14]];
    let d = t[15];
    let c = [t[3], t[7], t[11]];
    let e = Real::one() / (crate::vec3::dot(&c, &p) + d);
    let ee = e * e;
    let f = crate::vec3::add(&crate::mat3_col_major::mult_vec(&a, &p), &b);
    [
        a[0] * e - f[0] * c[0] * ee,
        a[1] * e - f[1] * c[0] * ee,
        a[2] * e - f[2] * c[0] * ee,
        a[3] * e - f[0] * c[1] * ee,
        a[4] * e - f[1] * c[1] * ee,
        a[5] * e - f[2] * c[1] * ee,
        a[6] * e - f[0] * c[2] * ee,
        a[7] * e - f[1] * c[2] * ee,
        a[8] * e - f[2] * c[2] * ee,
    ]
}

#[test]
fn test_jacobian_transform() {
    let a: [f64; 16] = [
        1.1, 2.3, 3.4, 1.4, 1.7, 3.2, -0.5, 0.2, 2.3, -1.3, 1.4, 0.3, 0.6, 2.3, 1.5, -2.3,
    ];
    let p0 = [1.3, 0.3, -0.5];
    let q0 = transform_homogeneous(&a, &p0).unwrap();
    let dqdp = jacobian_transform(&a, &p0);
    let eps = 1.0e-6;
    for i_dim in 0..3 {
        let p1 = {
            let mut p1 = p0;
            p1[i_dim] += eps;
            p1
        };
        let q1 = transform_homogeneous(&a, &p1).unwrap();
        for j_dim in 0..3 {
            let v_num = (q1[j_dim] - q0[j_dim]) / eps;
            let v_ana = dqdp[i_dim * 3 + j_dim];
            // dbg!(i_dim, j_dim, v_num, v_ana);
            assert!((v_num - v_ana).abs() < 9.0e-5);
        }
    }
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

pub fn try_inverse<Real>(b: &[Real; 16]) -> Option<[Real; 16]>
where
    Real: num_traits::Float + std::ops::MulAssign + std::ops::SubAssign,
{
    crate::matn::try_inverse::<Real, 4, 16>(b)
}

/// perspective transformation matrix (column major) compatible with blender
/// * asp - aspect ratio (width / height)
/// * lens - the focus distance (unit: mm) where the sensor size for longest edge is 18*2 mm.
/// * near - distance to the near clipping plane (>0)
/// * far - distance  to the far clipping plane (>0)
/// * proj_direction - projection direction
///     * true: projection from +Z
///     * false: projection from -Z
pub fn camera_perspective_blender<Real>(
    asp: Real,
    lens: Real,
    near: Real,
    far: Real,
    proj_direction: bool,
) -> [Real; 16]
where
    Real: num_traits::Float + 'static + std::ops::AddAssign,
    f64: AsPrimitive<Real>,
{
    if proj_direction {
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
    } else {
        let one = Real::one();
        let zero = Real::zero();
        let t0 = camera_perspective_blender(asp, lens, -far, -near, true);
        let cam_zflip: [Real; 16] = [
            -one, zero, zero, zero, zero, -one, zero, zero, zero, zero, -one, zero, zero, zero,
            zero, one,
        ];
        crate::mat4_col_major::mult_mat(&t0, &cam_zflip)
    }
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
    let transl = crate::mat4_col_major::from_translate(&[
        -cam_location[0],
        -cam_location[1],
        -cam_location[2],
    ]);
    let rot_x = crate::mat4_col_major::from_rot_x(-cam_rot_x_deg * deg2rad);
    let rot_y = crate::mat4_col_major::from_rot_y(-cam_rot_y_deg * deg2rad);
    let rot_z = crate::mat4_col_major::from_rot_z(-cam_rot_z_deg * deg2rad);
    let rot_yz = crate::mat4_col_major::mult_mat(&rot_y, &rot_z);
    let rot_zyx = crate::mat4_col_major::mult_mat(&rot_x, &rot_yz);
    crate::mat4_col_major::mult_mat(&rot_zyx, &transl)
}

pub fn scale<Real>(m: &[Real; 16], s: Real) -> [Real; 16]
where
    Real: Copy + std::ops::Mul<Output = Real>,
{
    m.map(|x| s * x)
}

pub fn mult_mat<Real>(a: &[Real; 16], b: &[Real; 16]) -> [Real; 16]
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
        1., 3., 4., 8., 3., 5., 5., 2., 5., 7., 8., 9., 8., 4., 5., 0.,
    ];
    let ainv = try_inverse(&a).unwrap();
    let ainv_a = mult_mat(&ainv, &a);
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

pub fn transpose<Real>(m: &[Real; 16]) -> [Real; 16]
where
    Real: Copy,
{
    [
        m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14], m[3], m[7],
        m[11], m[15],
    ]
}
