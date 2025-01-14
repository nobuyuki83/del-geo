//! symmetric matrix `[[a,b],[b,c]]` parameterized as `[a,b,c]`

pub fn mult_mat_sym<Real>(a: &[Real; 3], b: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    [
        a[0] * b[0] + a[1] * b[1],
        a[1] * b[0] + a[2] * b[1],
        a[0] * b[1] + a[1] * b[2],
        a[1] * b[1] + a[2] * b[2],
    ]
}

pub fn inverse<Real>(coeff: &[Real; 3]) -> Option<[Real; 3]>
where
    Real: num_traits::Float,
{
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    let det = a * c - b * b;
    if det.is_zero() {
        return None;
    }
    let di = Real::one() / det;
    Some([di * c, -di * b, di * a])
}

#[test]
fn test_inverse() {
    let a = [2.1, 0.1, 0.5];
    let ai = inverse(&a).unwrap();
    let aia = mult_mat_sym(&ai, &a);
    let diff = nalgebra::Vector4::<f32>::from_column_slice(&aia)
        - nalgebra::Vector4::<f32>::new(1., 0., 0., 1.);
    assert!(diff.norm() < 1.0e-7);
}

pub fn safe_inverse<Real>(coeff: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    let det = a * c - b * b;
    if det.abs() <= Real::epsilon() {
        let l = Real::one() / Real::epsilon();
        let a1 = a + Real::epsilon();
        let c1 = c + Real::epsilon();
        return [c1 * l, -b * l, a1 * l];
    }
    let di = Real::one() / det;
    [di * c, -di * b, di * a]
}

pub fn safe_inverse_preserve_positive_definiteness<Real>(abc: &[Real; 3], eps: Real) -> [Real; 3]
where
    Real: num_traits::Float + std::fmt::Display,
{
    assert!(
        abc[0] + abc[2] > Real::zero(),
        "{} {} {}",
        abc[0] + abc[2],
        abc[0],
        abc[2]
    );
    let eig_min = (abc[0] + abc[2]) * eps;
    if (abc[0] * abc[2] - abc[1] * abc[1]).abs() < eig_min {
        // one of the eigen value is nearly zero
        //let one = Real::one();
        let one = Real::one();
        let (e, v) = principal_directions(abc);
        // println!("　　　sig: {:?} {} {}",e, abc[0]*abc[2]-abc[1]*abc[1], abc[0]+abc[2]);
        let e0inv = one / (e[0] + eps);
        let e1inv = one / (e[1] + eps);
        [
            e0inv * v[0][0] * v[0][0] + e1inv * v[1][0] * v[1][0],
            e0inv * v[0][1] * v[0][0] + e1inv * v[1][1] * v[1][0],
            e0inv * v[0][1] * v[0][1] + e1inv * v[1][1] * v[1][1],
        ]
        // let (e,v) = del_geo_core::mat2_sym::prinsipal_directions(&xyz);
        // println!("　　　siginv: {:?}",e);
        //xyz
    } else {
        safe_inverse(abc)
    }
}

/// ax^2 + 2bxy + cy^2 = 1
pub fn principal_directions<Real>(coeff: &[Real; 3]) -> ([Real; 2], [[Real; 2]; 2])
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    if b.is_zero() {
        return ([a, c], [[one, zero], [zero, one]]);
    }
    let two = one + one;
    let four = two + two;
    let half = one / two;
    let tmp = ((a - c).powi(2) + four * b * b).sqrt();
    let lam0 = half * (a + c - tmp);
    let lam1 = half * (a + c + tmp);
    let det0 = a - c + tmp;
    let det1 = a - c - tmp;
    if det0.abs() > det1.abs() {
        let evec0 = [-two * b, det0];
        let evec1 = [det0, two * b];
        ([lam0, lam1], [evec0, evec1])
    } else {
        let evec0 = [det1, two * b];
        let evec1 = [-two * b, det1];
        ([lam0, lam1], [evec0, evec1])
    }
}

pub fn matvec<Real>(coeff: &[Real; 3], v: &[Real; 2]) -> [Real; 2]
where
    Real: num_traits::Float,
{
    [
        coeff[0] * v[0] + coeff[1] * v[1],
        coeff[1] * v[0] + coeff[2] * v[1],
    ]
}

pub fn aabb2<Real>(coeff: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let a = coeff[0];
    let b = coeff[1];
    let c = coeff[2];
    let det = a * c - b.powi(2);
    let minx = (c / det).sqrt();
    let miny = (a / det).sqrt();
    [-minx, -miny, minx, miny]
}

#[test]
fn test_prinsipal_directions() {
    let mats = [
        [4., 0., 1.],
        [4., 0.2, 1.],
        [1., 0., 4.],
        [1., -0.1, 4.],
        [4., 0., 4.],
        [4.0, 1.99, 1.0],
    ];
    for mat in mats {
        let ([lam0, lam1], [evec0, evec1]) = principal_directions::<f32>(&mat);
        let evec0 = nalgebra::Vector2::<f32>::from_row_slice(&evec0);
        let evec1 = nalgebra::Vector2::<f32>::from_row_slice(&evec1);
        assert!(evec0.norm() > 1.0e-5);
        assert!(evec1.norm() > 1.0e-5);
        let evec0 = evec0.normalize();
        let evec1 = evec1.normalize();
        assert!(evec0.dot(&evec1).abs() < 1.0e-5, "{}", evec0.dot(&evec1));
        {
            let tmp0 = nalgebra::Vector2::<f32>::from_row_slice(&matvec(
                &mat,
                evec0.as_slice()[..2].try_into().unwrap(),
            ));
            let diff = (tmp0 - lam0 * evec0).norm();
            assert!(diff < 1.0e-6, "{}", diff);
        }
        {
            let tmp0 = nalgebra::Vector2::<f32>::from_row_slice(&matvec(
                &mat,
                evec1.as_slice()[..2].try_into().unwrap(),
            ));
            assert!((tmp0 - lam1 * evec1).norm() < 1.0e-6);
        }
        let evec0 = evec0 / lam0.sqrt();
        let evec1 = evec1 / lam1.sqrt();
        let sig = nalgebra::Matrix2::<f32>::new(mat[0], mat[1], mat[1], mat[2]);
        let aabb = aabb2(&mat);
        let ndiv = 128;
        let mut sdf_max = f32::MIN;
        for i in 0..ndiv {
            let theta = (i as f32) * 2.0 * std::f32::consts::PI / (ndiv as f32);
            let v = evec0 * theta.cos() + evec1 * theta.sin();
            let radrad = (v.transpose() * sig * v)[0];
            assert!((radrad - 1.).abs() < 1.0e-3, "{}", radrad);
            assert!(crate::aabb2::is_include_point2(&aabb, &[v[0], v[1]]));
            let sdf = crate::aabb2::sdf(&aabb, &[v[0], v[1]]);
            assert!(sdf <= 0f32);
            sdf_max = sdf_max.max(sdf);
        }
        assert!(sdf_max > -1.0e-3, "{}", sdf_max);
    }
}

/// used for 3D Gaussian splatting
/// # return
/// sigma, Dsigma/Ddiarot (3x6)
pub fn wdw_projected_spd_mat3<Real>(
    p_mat: &[Real; 6],
    quat0: &[Real; 4],
    d: &[Real; 3],
) -> ([Real; 3], [[Real; 6]; 3])
where
    Real: num_traits::Float + std::ops::AddAssign,
{
    let two = Real::one() + Real::one();
    let r = crate::quaternion::to_mat3_col_major(quat0);
    let pr = crate::mat2x3_col_major::mult_mat3_col_major(p_mat, &r);
    let dd = [d[0] * d[0], d[1] * d[1], d[2] * d[2]];
    let a = pr[0] * pr[0] * dd[0] + pr[2] * pr[2] * dd[1] + pr[4] * pr[4] * dd[2];
    let b = pr[0] * pr[1] * dd[0] + pr[2] * pr[3] * dd[1] + pr[4] * pr[5] * dd[2];
    let c = pr[1] * pr[1] * dd[0] + pr[3] * pr[3] * dd[1] + pr[5] * pr[5] * dd[2];
    let mut dsdt = [[Real::zero(); 6]; 3];
    dsdt[0][0] = pr[0] * pr[0] * two * d[0];
    dsdt[1][0] = pr[0] * pr[1] * two * d[0];
    dsdt[2][0] = pr[1] * pr[1] * two * d[0];
    dsdt[0][1] = pr[2] * pr[2] * two * d[1];
    dsdt[1][1] = pr[2] * pr[3] * two * d[1];
    dsdt[2][1] = pr[3] * pr[3] * two * d[1];
    dsdt[0][2] = pr[4] * pr[4] * two * d[2];
    dsdt[1][2] = pr[4] * pr[5] * two * d[2];
    dsdt[2][2] = pr[5] * pr[5] * two * d[2];
    let rdd = [
        r[0] * dd[0],
        r[1] * dd[0],
        r[2] * dd[0],
        r[3] * dd[1],
        r[4] * dd[1],
        r[5] * dd[1],
        r[6] * dd[2],
        r[7] * dd[2],
        r[8] * dd[2],
    ];
    let rddrt = crate::mat3_col_major::mult_mat_row_major(&rdd, &r);
    let p0 = [p_mat[0], p_mat[2], p_mat[4]];
    let p1 = [p_mat[1], p_mat[3], p_mat[5]];
    let rddrtp0 = crate::mat3_col_major::mult_vec(&rddrt, &p0);
    let rddrtp1 = crate::mat3_col_major::mult_vec(&rddrt, &p1);
    let p0ddrtp0 = crate::vec3::cross(&p0, &rddrtp0);
    let p0ddrtp1 = crate::vec3::cross(&p0, &rddrtp1);
    let p1ddrtp0 = crate::vec3::cross(&p1, &rddrtp0);
    let p1ddrtp1 = crate::vec3::cross(&p1, &rddrtp1);
    for i in 0..3 {
        dsdt[0][3 + i] = -two * p0ddrtp0[i];
        dsdt[1][3 + i] = -p1ddrtp0[i] - p0ddrtp1[i];
        dsdt[2][3 + i] = -two * p1ddrtp1[i];
    }
    ([a, b, c], dsdt)
}

#[test]
fn test_wdw_projected_spd_mat3() {
    type Real = f64;
    let p_mat = [1., 2., 4., 3., 2., 0.];
    let quat0 = crate::quaternion::normalized(&[-3., -2., 0., -1.]);
    let s0_mat = [0.1, 3.0, 1.0];
    let (abc, dabcdt) = wdw_projected_spd_mat3(&p_mat, &quat0, &s0_mat);
    let p_mat = nalgebra::Matrix2x3::<Real>::from_column_slice(&p_mat);
    let r0_mat = crate::quaternion::to_mat3_col_major(&quat0);
    let r0_mat = nalgebra::Matrix3::<Real>::from_column_slice(&r0_mat);
    let s0_mat = crate::mat3_col_major::from_diagonal(&s0_mat);
    let s0_mat = nalgebra::Matrix3::<Real>::from_column_slice(&s0_mat);
    {
        let sigma0 = p_mat * r0_mat * s0_mat * s0_mat * r0_mat.transpose() * p_mat.transpose();
        assert!((abc[0] - sigma0.m11).abs() < 1.0e-5);
        assert!((abc[1] - sigma0.m12).abs() < 1.0e-5);
        assert!((abc[2] - sigma0.m22).abs() < 1.0e-5);
    }
    let eps: Real = 1.0e-5;
    for i in 0..3 {
        let mut s1_mat = s0_mat;
        s1_mat[(i, i)] += eps;
        let sigma1: nalgebra::Matrix2<Real> =
            p_mat * r0_mat * s1_mat * s1_mat.transpose() * r0_mat.transpose() * p_mat.transpose();
        let a_si_diff = (sigma1.m11 - abc[0]) / eps;
        let b_si_diff = (sigma1.m12 - abc[1]) / eps;
        let c_si_diff = (sigma1.m22 - abc[2]) / eps;
        let a_si_ana = dabcdt[0][i];
        let b_si_ana = dabcdt[1][i];
        let c_si_ana = dabcdt[2][i];
        assert!(
            (a_si_diff - a_si_ana).abs() < 2.0e-4,
            "{} {} {}",
            a_si_diff,
            a_si_ana,
            a_si_diff - a_si_ana
        );
        assert!(
            (b_si_diff - b_si_ana).abs() < 2.0e-4,
            "{} {} {}",
            b_si_diff,
            b_si_ana,
            b_si_diff - b_si_ana
        );
        assert!(
            (c_si_diff - c_si_ana).abs() < 2.0e-4,
            "{} {} {}",
            c_si_diff,
            c_si_ana,
            c_si_diff - c_si_ana
        );
    }
    for i in 0..3 {
        let aa = crate::vec3::basis(i, eps);
        let w = crate::vec3::to_mat3_from_axisangle_vec(&aa);
        let w = nalgebra::Matrix3::<Real>::from_column_slice(&w);
        let r1_mat = w * r0_mat;
        let sigma1 = p_mat * r1_mat * s0_mat * s0_mat * r1_mat.transpose() * p_mat.transpose();
        let a_ri_diff = (sigma1.m11 - abc[0]) / eps;
        let b_ri_diff = (sigma1.m12 - abc[1]) / eps;
        let c_ri_diff = (sigma1.m22 - abc[2]) / eps;
        let a_ri_ana = dabcdt[0][3 + i];
        let b_ri_ana = dabcdt[1][3 + i];
        let c_ri_ana = dabcdt[2][3 + i];
        dbg!(a_ri_diff, a_ri_ana);
        dbg!(c_ri_diff, c_ri_ana);
        dbg!(b_ri_diff, b_ri_ana);
    }
}

pub fn wdw_inverse<Real, const N: usize>(dabcdt: &[[Real; N]; 3], xyz: &[Real; 3]) -> [[Real; N]; 3]
where
    Real: num_traits::Float + std::ops::AddAssign + std::fmt::Debug,
{
    let one = Real::one();
    let two = one + one;
    let (x, y, z) = (xyz[0], xyz[1], xyz[2]);
    let mut dxyzdt = [[Real::zero(); N]; 3];
    for i in 0..N {
        let da = -dabcdt[0][i];
        let db = -dabcdt[1][i];
        let dc = -dabcdt[2][i];
        dxyzdt[0][i] = x * da * x + two * x * db * y + y * dc * y;
        dxyzdt[1][i] = x * da * y + z * db * x + y * db * y + y * dc * z;
        dxyzdt[2][i] = y * da * y + two * z * db * y + z * dc * z;
    }
    dxyzdt
}

#[test]
fn test_wdw_inv_projected_spd_mat3() {
    type Real = f64;
    let p_mat = [1., 2., 4., 3., 2., 0.];
    let quat0 = crate::quaternion::normalized(&[-3., -2., 0., -1.]);
    let s0_mat = [0.1, 3.0, 1.0];
    let (abc, dabcdt) = wdw_projected_spd_mat3(&p_mat, &quat0, &s0_mat);
    let xyz = safe_inverse(&abc);
    let dxyzdt = wdw_inverse(&dabcdt, &xyz);
    let p_mat = nalgebra::Matrix2x3::<Real>::from_column_slice(&p_mat);
    let r0_mat = crate::quaternion::to_mat3_col_major(&quat0);
    let r0_mat = nalgebra::Matrix3::<Real>::from_column_slice(&r0_mat);
    let s0_mat = crate::mat3_col_major::from_diagonal(&s0_mat);
    let s0_mat = nalgebra::Matrix3::<Real>::from_column_slice(&s0_mat);
    {
        let sigma0 = p_mat * r0_mat * s0_mat * s0_mat * r0_mat.transpose() * p_mat.transpose();
        let sigma0inv = sigma0.try_inverse().unwrap();
        assert!((xyz[0] - sigma0inv.m11).abs() < 1.0e-5);
        assert!((xyz[1] - sigma0inv.m12).abs() < 1.0e-5);
        assert!((xyz[2] - sigma0inv.m22).abs() < 1.0e-5);
    }
    let eps: Real = 1.0e-5;
    for i in 0..3 {
        let mut s1_mat = s0_mat;
        s1_mat[(i, i)] += eps;
        let sigma1: nalgebra::Matrix2<Real> =
            p_mat * r0_mat * s1_mat * s1_mat.transpose() * r0_mat.transpose() * p_mat.transpose();
        let sigma1inv = sigma1.try_inverse().unwrap();
        let x_si_diff = (sigma1inv.m11 - xyz[0]) / eps;
        let y_si_diff = (sigma1inv.m12 - xyz[1]) / eps;
        let z_si_diff = (sigma1inv.m22 - xyz[2]) / eps;
        let x_si_ana = dxyzdt[0][i];
        let y_si_ana = dxyzdt[1][i];
        let z_si_ana = dxyzdt[2][i];
        assert!(
            (x_si_diff - x_si_ana).abs() < 1.0e-3,
            "{} {} {}",
            x_si_diff,
            x_si_ana,
            x_si_diff - x_si_ana
        );
        assert!(
            (y_si_diff - y_si_ana).abs() < 1.0e-3,
            "{} {} {}",
            y_si_diff,
            y_si_ana,
            y_si_diff - y_si_ana
        );
        assert!(
            (z_si_diff - z_si_ana).abs() < 1.0e-3,
            "{} {} {}",
            z_si_diff,
            z_si_ana,
            z_si_diff - z_si_ana
        );
    }
    for i in 0..3 {
        let aa = crate::vec3::basis(i, eps);
        let w = crate::vec3::to_mat3_from_axisangle_vec(&aa);
        let w = nalgebra::Matrix3::<Real>::from_column_slice(&w);
        let r1_mat = w * r0_mat;
        let sigma1 = p_mat * r1_mat * s0_mat * s0_mat * r1_mat.transpose() * p_mat.transpose();
        let sigma1inv = sigma1.try_inverse().unwrap();
        let x_ri_diff = (sigma1inv.m11 - xyz[0]) / eps;
        let y_ri_diff = (sigma1inv.m12 - xyz[1]) / eps;
        let z_ri_diff = (sigma1inv.m22 - xyz[2]) / eps;
        let x_ri_ana = dxyzdt[0][3 + i];
        let y_ri_ana = dxyzdt[1][3 + i];
        let z_ri_ana = dxyzdt[2][3 + i];
        assert!(
            (x_ri_diff - x_ri_ana).abs() < 1.0e-3,
            "{} {} {}",
            x_ri_diff,
            x_ri_ana,
            x_ri_diff - x_ri_ana
        );
        assert!(
            (y_ri_diff - y_ri_ana).abs() < 1.0e-3,
            "{} {} {}",
            y_ri_diff,
            y_ri_ana,
            y_ri_diff - y_ri_ana
        );
        assert!(
            (z_ri_diff - z_ri_ana).abs() < 1.0e-3,
            "{} {} {}",
            z_ri_diff,
            z_ri_ana,
            z_ri_diff - z_ri_ana
        );
    }
}

pub fn mult_vec_from_both_sides<Real>(m: &[Real; 3], b: &[Real; 2], c: &[Real; 2]) -> Real
where
    Real: num_traits::Float,
{
    m[0] * b[0] * c[0] + m[1] * (b[0] * c[1] + b[1] * c[0]) + m[2] * b[1] * c[1]
}
