//! the trait and methods for 4x4 matrix with column major storage

/// trait for 4x4 matrix
pub trait Mat4ColMajor<Real>
where
    Self: Sized,
{
    /// transform a 3D point as homogeneous coordinate `(x,y,z,1)` and divide by `w`.
    /// return `(transformed point, w)`, or `None` if `w` is zero
    fn transform_homogeneous(&self, a: &[Real; 3]) -> Option<([Real; 3], Real)>;
    /// matrix product `self * b`
    fn mult_mat(&self, b: &Self) -> Self;
    /// transform a direction vector using only the upper-left 3x3 block (translation ignored)
    fn transform_direction(&self, a: &[Real; 3]) -> [Real; 3];
    /// inverse matrix, or `None` if the matrix is singular
    fn try_inverse(&self) -> Option<Self>;
    /// element-wise addition `self += other`
    fn add_in_place(&mut self, other: &Self);
    /// transposed matrix
    fn transpose(&self) -> [Real; 16];
}

impl<Real> Mat4ColMajor<Real> for [Real; 16]
where
    Real: num_traits::Float,
{
    fn transform_homogeneous(&self, v: &[Real; 3]) -> Option<([Real; 3], Real)> {
        transform_homogeneous(self, v)
    }
    fn mult_mat(&self, b: &Self) -> Self {
        mult_mat_col_major(self, b)
    }
    fn transform_direction(&self, a: &[Real; 3]) -> [Real; 3] {
        transform_direction(self, a)
    }
    fn try_inverse(&self) -> Option<Self> {
        try_inverse(self)
    }
    fn add_in_place(&mut self, other: &Self) {
        self.iter_mut()
            .zip(other.iter())
            .for_each(|(a, &b)| *a = *a + b);
    }
    fn transpose(&self) -> [Real; 16] {
        transpose(self)
    }
}

use crate::aabb3::max_edge_size;
use num_traits::AsPrimitive;

/// 4x4 identity matrix
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

/// diagonal matrix `diag(m11, m22, m33, m44)`
pub fn from_diagonal<Real>(m11: Real, m22: Real, m33: Real, m44: Real) -> [Real; 16]
where
    Real: num_traits::Zero + Copy,
{
    let zero = Real::zero();
    [
        m11, zero, zero, zero, zero, m22, zero, zero, zero, zero, m33, zero, zero, zero, zero, m44,
    ]
}

/// homogeneous transformation of uniform scaling by `s` (i.e., `diag(s, s, s, 1)`)
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

/// homogeneous transformation of translation by `v`
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

/// homogeneous transformation of rotation around the x-axis by `theta` (radian, right-handed)
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

/// homogeneous transformation of rotation around the y-axis by `theta` (radian, right-handed)
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

/// homogeneous transformation of rotation around the z-axis by `theta` (radian, right-handed)
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

/// rotation matrix where x-rotation, y-rotation and z-rotation is applied sequentially
/// (i.e., `Rz(rz) * Ry(ry) * Rx(rx)`). angles are in radian
pub fn from_bryant_angles<Real>(rx: Real, ry: Real, rz: Real) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let x = from_rot_x(rx);
    let y = from_rot_y(ry);
    let z = from_rot_z(rz);
    let yx = mult_mat_col_major(&y, &x);
    mult_mat_col_major(&z, &yx)
}

/// transformation converting normalized device coordinate (NDC) `[-1,+1]^3` to pixel coordinate
/// depth (-1, +1) is transformed to (0, +1)
/// for example:
///     [-1,-1,-1] becomes (0, H, 0)
///     [+1,+1,+1] becomes (W, 0, 1)
///
/// * Arguments
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

/*
/// transformation that maps the NDC cube `[-1,+1]^3` to a box centered at the AABB that covers it,
/// where the x/y aspect ratio of the box is `asp` (width / height).
/// The inverse of this matrix fits the AABB into the NDC while preserving the x/y aspect ratio.
/// * `aabb` - `[min_x, min_y, min_z, max_x, max_y, max_z]`
/// * `asp` - aspect ratio (width / height) of the screen
pub fn from_aabb3_fit_into_ndc_preserving_xyasp<Real>(aabb: &[Real; 6], asp: Real) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    let cntr = crate::aabb3::center(aabb);
    let (scale_xy, scale_z) = {
        let size = crate::aabb3::size(aabb);
        let lenx = size[0] / asp;
        let leny = size[1];
        if lenx > leny {
            (lenx / two, size[2] / two)
        } else {
            (leny / two, size[2] / two)
        }
    };
    [
        scale_xy * asp,
        zero,
        zero,
        zero,
        zero,
        scale_xy,
        zero,
        zero,
        zero,
        zero,
        scale_z,
        zero,
        cntr[0],
        cntr[1],
        cntr[2],
        one,
    ]
}
 */

/// transform aabb to unit square (0,1)^3 while preserving aspect ratio
/// (uniformly scaled so that the longest edge becomes 1, and centered at (0.5, 0.5, 0.5))
/// return 4x4 homogeneous transformation matrix in **column major** order
pub fn from_aabb3_fit_into_unit_preserve_asp<Real>(aabb_world: &[Real; 6]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    let half = one / two;
    let cntr = [
        (aabb_world[0] + aabb_world[3]) / two,
        (aabb_world[1] + aabb_world[4]) / two,
        (aabb_world[2] + aabb_world[5]) / two,
    ];
    let size = max_edge_size(aabb_world);
    let a = one / size;
    let b = -a * cntr[0] + half;
    let c = -a * cntr[1] + half;
    let d = -a * cntr[2] + half;
    [
        a, zero, zero, zero, zero, a, zero, zero, zero, zero, a, zero, b, c, d, one,
    ]
}

/// transform aabb to unit square (0,1)^3 (each axis is scaled independently)
/// return 4x4 homogeneous transformation matrix in **column major** order
pub fn from_aabb3_fit_into_unit<Real>(aabb_world: &[Real; 6]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    let half = one / two;
    let cntr = crate::aabb3::center(aabb_world);
    let size = crate::aabb3::size(aabb_world);
    let ax = one / size[0];
    let ay = one / size[1];
    let az = one / size[2];
    let b = -ax * cntr[0] + half;
    let c = -ay * cntr[1] + half;
    let d = -az * cntr[2] + half;
    [
        ax, zero, zero, zero, zero, ay, zero, zero, zero, zero, az, zero, b, c, d, one,
    ]
}

/// this function is typically used to make 3D homogeneous tranformation matrix
/// from 2D homogeneous transformation mtrix.
/// x and y are transformed by `m` while z is left unchanged
pub fn from_mat3_col_major_adding_z<Real>(m: &[Real; 9]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    [
        m[0], m[1], zero, m[2], m[3], m[4], zero, m[5], zero, zero, one, zero, m[6], m[7], zero,
        m[8],
    ]
}

/// embed a 3x3 matrix `m` to the upper-left block of a 4x4 matrix.
/// the bottom-right element is set to `dia`, and the other elements are zero
pub fn from_mat3_col_major_adding_w<Real>(m: &[Real; 9], dia: Real) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    [
        m[0], m[1], m[2], zero, m[3], m[4], m[5], zero, m[6], m[7], m[8], zero, zero, zero, zero,
        dia,
    ]
}

// above: from method (making 4x4 matrix)
// ----------------------------------------

/// extract the upper-left 3x3 block (column major)
pub fn to_mat3_col_major_xyz<T>(m: &[T; 16]) -> [T; 9]
where
    T: num_traits::Float,
{
    [m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9], m[10]]
}

/// extract the translation part (the first three elements of the fourth column)
pub fn to_vec3_translation<T>(m: &[T; 16]) -> [T; 3]
where
    T: num_traits::Float,
{
    [m[12], m[13], m[14]]
}

// above: to method
// ----------------------------------------

/// transform a 3D point `x` as homogeneous coordinate `(x,y,z,1)` and divide by `w`.
/// return `(transformed point, w)`, or `None` if `w` is zero
pub fn transform_homogeneous<Real>(
    transform: &[Real; 16],
    x: &[Real; 3],
) -> Option<([Real; 3], Real)>
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
    Some(([y0 / y3, y1 / y3, y2 / y3], y3))
}

/// jacobian of the homogeneous transformation (see `transform_homogeneous`) at point `p`.
/// return 3x3 matrix in column major order where `(i,j)` element is `d q_i / d p_j`
pub fn jacobian_transform<Real>(t: &[Real; 16], p: &[Real; 3]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    let a = [t[0], t[1], t[2], t[4], t[5], t[6], t[8], t[9], t[10]];
    let b = [t[12], t[13], t[14]];
    let d = t[15];
    let c = [t[3], t[7], t[11]];
    let e = Real::one() / (crate::vec3::dot(&c, p) + d);
    let ee = e * e;
    let f = crate::vec3::add(&crate::mat3_col_major::mult_vec(&a, p), &b);
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
    let a0: [f64; 16] = [
        1.1, 2.3, 3.4, 1.4, 1.7, 3.2, -0.5, 0.2, 2.3, -1.3, 1.4, 0.3, 0.6, 2.3, 1.5, -2.3,
    ];
    let a1 = [
        -0.9, 0.0, 0.0, 0.0, 0.0, -0.9, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.5, -2.0,
    ];
    let vec_a = [a0, a1];
    for a in vec_a.iter() {
        // let p0 = [1.3, 0.3, -0.5];
        let p0 = [0.5, 0.3, 0.0];
        let (q0, _) = transform_homogeneous(a, &p0).unwrap();
        let dqdp = jacobian_transform(a, &p0);
        let eps = 1.0e-6;
        for j_dim in 0..3 {
            let p1 = {
                let mut p1 = p0;
                p1[j_dim] += eps;
                p1
            };
            let (q1, _) = transform_homogeneous(a, &p1).unwrap();
            for i_dim in 0..3 {
                let v_num = (q1[i_dim] - q0[i_dim]) / eps;
                let v_ana = dqdp[i_dim + 3 * j_dim];
                // dbg!(i_dim, j_dim, v_num, v_ana);
                assert!((v_num - v_ana).abs() < 9.0e-5);
            }
        }
    }
}

/// transform a direction vector `x` using only the upper-left 3x3 block (translation ignored)
pub fn transform_direction<Real>(transform: &[Real; 16], x: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let y0 = transform[0] * x[0] + transform[4] * x[1] + transform[8] * x[2];
    let y1 = transform[1] * x[0] + transform[5] * x[1] + transform[9] * x[2];
    let y2 = transform[2] * x[0] + transform[6] * x[1] + transform[10] * x[2];
    [y0, y1, y2]
}

/// inverse matrix, or `None` if the matrix is singular.
/// (the row-major routine can be used because `inv(A^T) = inv(A)^T`)
pub fn try_inverse<Real>(b: &[Real; 16]) -> Option<[Real; 16]>
where
    Real: num_traits::Float,
{
    crate::matn_row_major::try_inverse::<Real, 4, 16>(b)
}

/// inverse matrix computed by LU decomposition with partial pivoting.
/// more robust than `try_inverse` when a diagonal element is (near) zero.
/// return `None` if a pivot is smaller than the machine epsilon
pub fn try_inverse_with_pivot<Real>(a: &[Real; 16]) -> Option<[Real; 16]>
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    // --- Step 1: Convert to row-major scratch array for LU ---
    let mut piv = [0usize; 4];
    for i in 0..4 {
        piv[i] = i;
    }
    let eps = Real::epsilon();
    let mut a = a.to_owned();

    // --- LU decomposition with partial pivoting ---
    for k in 0..4 {
        // pivot row
        let mut p = k;
        let mut maxv = a[k * 4 + k].abs();
        for r in (k + 1)..4 {
            let v = a[r * 4 + k].abs();
            if v > maxv {
                maxv = v;
                p = r;
            }
        }
        if maxv <= eps {
            return None;
        }
        if p != k {
            for c in 0..4 {
                a.swap(k * 4 + c, p * 4 + c);
            }
            piv.swap(k, p);
        }
        let pivot = a[k * 4 + k];
        for r in (k + 1)..4 {
            a[r * 4 + k] = a[r * 4 + k] / pivot;
            let m = a[r * 4 + k];
            for c in (k + 1)..4 {
                a[r * 4 + c] = a[r * 4 + c] - m * a[k * 4 + c];
            }
        }
    }

    // --- Helper: solve (LU)x = P*b ---
    let solve = |b: [Real; 4]| -> [Real; 4] {
        let mut y = [zero; 4];
        for i in 0..4 {
            y[i] = b[piv[i]];
        }
        // forward
        for i in 0..4 {
            for j in 0..i {
                y[i] = y[i] - a[i * 4 + j] * y[j];
            }
        }
        // back
        let mut x = [zero; 4];
        for i in (0..4).rev() {
            let mut s = y[i];
            for j in (i + 1)..4 {
                s = s - a[i * 4 + j] * x[j];
            }
            x[i] = s / a[i * 4 + i];
        }
        x
    };

    // --- Step 2: Solve for each column of identity ---
    let mut inv_cols = [zero; 16];
    for col in 0..4 {
        let mut e = [zero; 4];
        e[col] = one;
        let x = solve(e);
        for row in 0..4 {
            inv_cols[col + 4 * row] = x[row];
        }
    }
    Some(inv_cols)
}

#[test]
fn test_try_inverse_with_pivot() {
    let a = [
        1.0, 3.0, 2.0, 4.0, 2.0, 5.0, 3.0, 1.0, 0.0, 1.0, 4.0, 3.0, 4.0, 2.0, 1.0, 3.0,
    ];
    let ainv = try_inverse_with_pivot(&a).unwrap();
    let a_ainv = mult_mat_col_major(&a, &ainv);
    dbg!(a_ainv);
    //
    let a = [
        0., 0., 1., 0., -1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.,
    ];
    let ainv = try_inverse_with_pivot(&a).unwrap();
    let a_ainv = mult_mat_col_major(&a, &ainv);
    dbg!(a_ainv);
}

/// perspective transformation matrix (column major) compatible with blender
/// * asp - aspect ratio (width / height)
/// * lens - the focus distance (unit: mm) where the sensor size for longest edge is 18*2 mm.
/// * near - distance to the near clipping plane (>0)
/// * far - distance  to the far clipping plane (>0)
/// * proj_direction - projection direction
///     * true: projection from +Z (ray is going -Z direction, compatible to this library's definition)
///     * false: projection from -Z (ray is going +Z direction, compatible to OpenGl)
pub fn camera_perspective_blender<Real>(
    asp: Real,
    lens: Real,
    near: Real,
    far: Real,
    proj_direction: bool,
) -> [Real; 16]
where
    Real: num_traits::Float + 'static,
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
        crate::mat4_col_major::mult_mat_col_major(&t0, &cam_zflip)
    }
}

/// view matrix (transformation from world to camera coordinate) compatible with blender
/// * `cam_location` - location of the camera in the world coordinate
/// * `cam_rot_{x,y,z}_deg` - XYZ euler angles of the camera (unit: degree)
pub fn camera_external_blender<Real>(
    cam_location: &[Real; 3],
    cam_rot_x_deg: Real,
    cam_rot_y_deg: Real,
    cam_rot_z_deg: Real,
) -> [Real; 16]
where
    Real: num_traits::Float + num_traits::FloatConst + 'static,
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
    let rot_yz = crate::mat4_col_major::mult_mat_col_major(&rot_y, &rot_z);
    let rot_zyx = crate::mat4_col_major::mult_mat_col_major(&rot_x, &rot_yz);
    crate::mat4_col_major::mult_mat_col_major(&rot_zyx, &transl)
}

/// multiply all the elements by the scalar `s`
pub fn scale<Real>(m: &[Real; 16], s: Real) -> [Real; 16]
where
    Real: Copy + std::ops::Mul<Output = Real>,
{
    m.map(|x| s * x)
}

/// matrix product `a * b`
pub fn mult_mat_col_major<Real>(a: &[Real; 16], b: &[Real; 16]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let mut o = [Real::zero(); 16];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                o[i + j * 4] = o[i + j * 4] + a[i + k * 4] * b[k + j * 4];
            }
        }
    }
    o
}

/// matrix-vector product `a * b` for a 4D vector `b`
pub fn mult_vec<Real>(a: &[Real; 16], b: &[Real; 4]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let mut o = [Real::zero(); 4];
    for i in 0..4 {
        for j in 0..4 {
            o[i] = o[i] + a[i + 4 * j] * b[j];
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
    let ainv_a = mult_mat_col_major(&ainv, &a);
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

/// transposed matrix
pub fn transpose<Real>(m: &[Real; 16]) -> [Real; 16]
where
    Real: Copy,
{
    [
        m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14], m[3], m[7],
        m[11], m[15],
    ]
}

/// ray that goes through `pos_world: [f32;3]` that will be -z direction in the normalized device coordinate (NDC).
/// the ray starts on the front plane (z=+1 in NDC) and ends on the back plane (z=-1 in NDC).
/// return `(ray_org: [f32;3], ray_dir: [f32;3])`
pub fn ray_from_transform_world2ndc(
    transform_world2ndc: &[f32; 16],
    pos_world: &[f32; 3],
    transform_ndc2world: &[f32; 16],
) -> ([f32; 3], [f32; 3]) {
    let (pos_mid_ndc, _) = transform_homogeneous(transform_world2ndc, pos_world).unwrap();
    let (ray_stt_world, _) =
        transform_homogeneous(transform_ndc2world, &[pos_mid_ndc[0], pos_mid_ndc[1], 1.0]).unwrap();
    let (ray_end_world, _) =
        transform_homogeneous(transform_ndc2world, &[pos_mid_ndc[0], pos_mid_ndc[1], -1.0])
            .unwrap();
    (
        ray_stt_world,
        crate::vec3::sub(&ray_end_world, &ray_stt_world),
    )
}

/// the ray start from the front plane and ends on the back plane
/// "integer corner" coordinate is used here (the center of the top-left pixel is (0.5, 0.5) )
pub fn ray_from_transform_ndc2world_and_pixel_coordinates<Real>(
    pix_coord: (Real, Real),
    image_size: &(Real, Real),
    transform_ndc2world: &[Real; 16],
) -> ([Real; 3], [Real; 3])
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let x0 = two * pix_coord.0 / (image_size.0) - one;
    let y0 = one - two * pix_coord.1 / (image_size.1);
    let (p0, _) = transform_homogeneous(transform_ndc2world, &[x0, y0, one]).unwrap();
    let (p1, _) = transform_homogeneous(transform_ndc2world, &[x0, y0, -one]).unwrap();
    let ray_org = p0;
    let ray_dir = crate::vec3::sub(&p1, &p0);
    (ray_org, ray_dir)
}

/// matrix product of three matrices `a * b * c`
pub fn mult_three_mats_col_major<Real>(a: &[Real; 16], b: &[Real; 16], c: &[Real; 16]) -> [Real; 16]
where
    Real: num_traits::Float,
{
    let d = mult_mat_col_major(b, c);
    mult_mat_col_major(a, &d)
}
