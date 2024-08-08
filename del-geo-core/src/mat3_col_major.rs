//! methods for 3x3 matrix where storage is column major order

use std::ops::AddAssign;

pub fn from_diagonal<Real>(s: &[Real; 3]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    [s[0], zero, zero, zero, s[1], zero, zero, zero, s[2]]
}

pub fn from_column_vectors<Real>(x: &[Real;3], y: &[Real;3], z: &[Real;3]) -> [Real;9]
where Real: Copy
{
    [
        x[0], x[1], x[2],
        y[0], y[1], y[2],
        z[0], z[1], z[2]
    ]
}

pub fn try_inverse<T>(b: &[T; 9]) -> Option<[T; 9]>
where
    T: num_traits::Float,
{
    let det = b[0] * b[4] * b[8] + b[3] * b[7] * b[2] + b[6] * b[1] * b[5]
        - b[0] * b[7] * b[5]
        - b[6] * b[4] * b[2]
        - b[3] * b[1] * b[8];
    if det.is_zero() {
        return None;
    }
    let inv_det = T::one() / det;
    Some([
        inv_det * (b[4] * b[8] - b[5] * b[7]),
        inv_det * (b[2] * b[7] - b[1] * b[8]),
        inv_det * (b[1] * b[5] - b[2] * b[4]),
        inv_det * (b[5] * b[6] - b[3] * b[8]),
        inv_det * (b[0] * b[8] - b[2] * b[6]),
        inv_det * (b[2] * b[3] - b[0] * b[5]),
        inv_det * (b[3] * b[7] - b[4] * b[6]),
        inv_det * (b[1] * b[6] - b[0] * b[7]),
        inv_det * (b[0] * b[4] - b[1] * b[3]),
    ])
}

pub fn transform_homogeneous<Real>(transform: &[Real; 9], x: &[Real; 2]) -> Option<[Real; 2]>
where
    Real: num_traits::Float,
{
    let y2 = transform[2] * x[0] + transform[5] * x[1] + transform[8];
    if y2.is_zero() {
        return None;
    }
    //
    let y0 = transform[0] * x[0] + transform[3] * x[1] + transform[6];
    let y1 = transform[1] * x[0] + transform[4] * x[1] + transform[7];
    Some([y0 / y2, y1 / y2])
}

pub fn identity<T>() -> [T; 9]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, zero, one, zero, zero, zero, one]
}

pub fn mult_vec<Real>(a: &[Real; 9], b: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    [
        a[0] * b[0] + a[3] * b[1] + a[6] * b[2],
        a[1] * b[0] + a[4] * b[1] + a[7] * b[2],
        a[2] * b[0] + a[5] * b[1] + a[8] * b[2],
    ]
}

pub fn transpose<Real>(m: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

pub fn mult_mat_col_major<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float + AddAssign,
{
    let mut r = [Real::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i + 3 * j] += a[i + 3 * k] * b[k + 3 * j]
            }
        }
    }
    r
}

pub fn mult_mat_row_major<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
    where
        Real: num_traits::Float + AddAssign,
{
    let mut r = [Real::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i + 3 * j] += a[i + 3 * k] * b[j + 3 * k]
            }
        }
    }
    r
}

pub fn to_quaternion<Real>(p: &[Real; 9]) -> [Real; 4]
where
    Real: num_traits::Float + ordered_float::FloatCore + std::fmt::Debug,
{
    let one = Real::one();
    let one4th = one / (one + one + one + one);
    let smat: [Real; 16] = [
        one + p[0] - p[4] - p[8], // 00
        p[3] + p[1],              // 01
        p[6] + p[2],              // 02
        p[5] - p[7],              // 03
        p[1] + p[3],              // 10
        one - p[0] + p[4] - p[8], // 11
        p[7] + p[5],              // 12
        p[6] - p[2],              // 13
        p[6] + p[2],              // 20
        p[7] + p[5],              // 21
        one - p[0] - p[4] + p[8], // 22
        p[1] - p[3],              // 23
        p[5] - p[7],              // 30
        p[6] - p[2],              // 31
        p[1] - p[3],              // 32
        one + p[0] + p[4] + p[8], // 33
    ];

    let dias = [smat[0], smat[5], smat[10], smat[15]];
    let imax = (0usize..4usize)
        .max_by_key(|&a| ordered_float::OrderedFloat(dias[a]))
        .unwrap();
    assert!(dias[0] <= dias[imax]);
    assert!(dias[1] <= dias[imax]);
    assert!(dias[2] <= dias[imax]);
    assert!(dias[3] <= dias[imax]);

    let mut quat = [Real::zero(); 4];
    quat[imax] = smat[imax * 4 + imax].sqrt() / (one + one);
    for k in 0..4 {
        if k == imax {
            continue;
        } else {
            quat[k] = smat[imax * 4 + k] * one4th / quat[imax];
        }
    }
    return quat;
}

#[test]
fn test_to_quaternion() {
    let quats = [
        [-3., -2., 0., -1.],
        [3., -2., 0., -1.],
        [-1., 3., -2., -1.],
        [-1., -3., -2., -1.],
        [-1., -2., 3., -1.],
        [-1., -2., -3., -1.],
        [-1., -2., 1., -4.],
        [-1., -2., -1., -4.],
    ];
    for quat0 in quats {
        let quat0 = crate::quat::normalized(&quat0);
        let r_mat = crate::quat::to_mat3_col_major(&quat0);
        let quat1 = to_quaternion(&r_mat);
        let quat0 = nalgebra::Vector4::<f32>::from_row_slice(&quat0);
        let quat1 = nalgebra::Vector4::<f32>::from_row_slice(&quat1);
        assert!((quat0 - quat1).norm().min((quat0 + quat1).norm()) < 1.0e-7);
    }
}

// https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
pub fn to_vec3_axisangle_from_rot_mat<T>(m: [T; 9]) -> [T; 3]
where
    T: num_traits::Float,
{
    let one = T::one();
    let half = one / (one + one);
    let cos_t0 = (m[0] + m[4] + m[8] - one) * half;
    if (cos_t0 - one).abs() <= T::epsilon() {
        // very small rotation
        return [
            (m[5] - m[7]) * half,
            (m[6] - m[2]) * half,
            (m[1] - m[3]) * half,
        ];
    }
    let t0 = cos_t0.acos();
    let c0 = t0 * half / t0.sin();
    [c0 * (m[5] - m[7]), c0 * (m[6] - m[2]), c0 * (m[1] - m[3])]
}
