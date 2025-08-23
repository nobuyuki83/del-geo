//! symmetric matrix `[[a,b],[b,c]]` parameterized as `[a,b,c]`

/// symmetric matrix `[[a,b],[b,c]]` parameterized as `[a,b,c]`
pub trait Mat2Sym<Real>
where
    Self: Sized,
{
    fn mult_mat_sym(&self, other: &Self) -> [Real; 4];
    fn inverse(&self) -> Option<Self>;
    fn safe_inverse(&self) -> Self;
    fn eigen_decomposition(&self) -> ([Real; 4], [Real; 2]);
    fn mult_vec(&self, v: &[Real; 2]) -> [Real; 2];
    fn aabb2(&self) -> [Real; 4];
}

impl<Real> Mat2Sym<Real> for [Real; 3]
where
    Real: num_traits::Float + std::fmt::Debug,
{
    fn mult_mat_sym(&self, other: &Self) -> [Real; 4] {
        mult_mat_sym(self, other)
    }
    fn inverse(&self) -> Option<Self> {
        inverse(self)
    }
    fn safe_inverse(&self) -> Self {
        safe_inverse(self)
    }
    fn eigen_decomposition(&self) -> ([Real; 4], [Real; 2]) {
        eigen_decomposition(self)
    }
    fn mult_vec(&self, v: &[Real; 2]) -> [Real; 2] {
        mult_vec(self, v)
    }
    fn aabb2(&self) -> [Real; 4] {
        aabb2(self)
    }
}

// ----------------------------
// below: to methods

pub fn to_mat_col_major<Real>(sm: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    [sm[0], sm[1], sm[1], sm[2]]
}

// above: to methods
// --------------------------

pub fn from_columns<Real>(x: &[Real; 2], y: &[Real; 2]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    [x[0], x[1], y[0], y[1]]
}

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
    let ai = a.inverse().unwrap();
    let aia = ai.mult_mat_sym(&a);
    let diff = crate::vec4::distance(&aia, &[1., 0., 0., 1.]);
    assert!(diff < 1.0e-7);
}

pub fn safe_inverse<Real>(&[a, b, c]: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
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

pub fn safe_inverse_preserve_positive_definiteness<Real>(
    &[a, b, c]: &[Real; 3],
    eps: Real,
) -> [Real; 3]
where
    Real: num_traits::Float + std::fmt::Display + std::fmt::Debug,
{
    assert!(a + c > Real::zero(), "{} {} {}", a + c, a, c);
    let eig_min = (a + c) * eps;
    if (a * c - b * b).abs() < eig_min {
        // one of the eigen value is nearly zero
        //let one = Real::one();
        let one = Real::one();
        let (v, e) = [a, b, c].eigen_decomposition();
        let (v0, v1) = crate::mat2_col_major::to_columns(&v);
        // println!("　　　sig: {:?} {} {}",e, a*c-b*b, a+c);
        let e0inv = one / (e[0] + eps);
        let e1inv = one / (e[1] + eps);
        [
            e0inv * v0[0] * v0[0] + e1inv * v1[0] * v1[0],
            e0inv * v0[1] * v0[0] + e1inv * v1[1] * v1[0],
            e0inv * v0[1] * v0[1] + e1inv * v1[1] * v1[1],
        ]
        // let (e,v) = del_geo_core::mat2_sym::prinsipal_directions(&xyz);
        // println!("　　　siginv: {:?}",e);
        //xyz
    } else {
        [a, b, c].safe_inverse()
    }
}

/// this function returns U and \Sigma
/// A = U * \Sigma * U^t
/// the axis of ellipsoid: ax^2 + 2bxy + cy^2 = 1
pub fn eigen_decomposition<Real>(&[a, b, c]: &[Real; 3]) -> ([Real; 4], [Real; 2])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let zero = Real::zero();
    let one = Real::one();
    if b.is_zero() {
        return ([one, zero, zero, one], [a, c]);
    }
    let two = one + one;
    let four = two + two;
    let half = one / two;
    let tmp = ((a - c).powi(2) + four * b * b).sqrt();
    let lam0 = half * (a + c - tmp);
    let lam1 = half * (a + c + tmp);
    let det0 = a - c + tmp;
    let det1 = a - c - tmp;
    let u0 = if det0.abs() > det1.abs() {
        [-two * b, det0]
    } else {
        [det1, two * b]
    };
    let u0 = crate::vec2::normalize(&u0);
    let u1 = crate::vec2::rotate90(&u0);
    let u = from_columns(&u0, &u1);
    (u, [lam0, lam1])
}

#[test]
fn test_eigen_decomposition() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    // std::uniform_real_distribution < double > dist(-50.0, 50.0);
    for _itr in 0..1000 {
        let sm = std::array::from_fn(|_| rng.random_range(-30f64..30f64));
        let (u, l) = eigen_decomposition(&sm);
        use crate::mat2_col_major::Mat2ColMajor;
        let ut = u.transpose();
        {
            let utu = u.mult_mat_col_major(&ut);
            let err = crate::mat2_col_major::from_identity()
                .sub(&utu)
                .squared_norm();
            assert!(err < 1.0e-20);
        }
        {
            // check determinant of `u` is one
            let det = u.determinant();
            assert!((det - 1.).abs() < 1.0e-10, "{}", (det - 1.0).abs());
        }
        let l = crate::mat2_col_major::from_diagonal(&l);
        let ulut = u.mult_mat_col_major(&l).mult_mat_col_major(&ut);
        let sm = to_mat_col_major(&sm);
        let err = sm.sub(&ulut).squared_norm();
        assert!(err < 1.0e-10);
    }
}

pub fn mult_vec<Real>(&[c0, c1, c2]: &[Real; 3], &[v0, v1]: &[Real; 2]) -> [Real; 2]
where
    Real: num_traits::Float,
{
    [c0 * v0 + c1 * v1, c1 * v0 + c2 * v1]
}

/// 2-dimensional AABB for the 2D ellipsoid
/// (x,y)^T S (x,y) = 1
pub fn aabb2<Real>(&[a, b, c]: &[Real; 3]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let det = a * c - b.powi(2);
    let minx = (c / det).sqrt();
    let miny = (a / det).sqrt();
    [-minx, -miny, minx, miny]
}

pub fn determinant<Real>(&[a, b, c]: &[Real; 3]) -> Real
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    a * c - two * b
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
        use crate::vec2::Vec2;
        let (u_mat, [lam0, lam1]) = eigen_decomposition::<f32>(&mat);
        let (evec0, evec1) = crate::mat2_col_major::to_columns(&u_mat);
        assert!(crate::vec2::length(&evec0) > 1.0e-5);
        assert!(crate::vec2::length(&evec1) > 1.0e-5);
        let evec0 = evec0.normalize();
        let evec1 = evec1.normalize();
        assert!(evec0.dot(&evec1).abs() < 1.0e-5, "{}", evec0.dot(&evec1));
        {
            let tmp0 = &mult_vec(&mat, evec0.as_slice()[..2].try_into().unwrap());
            assert!(tmp0.sub(&evec0.scale(lam0)).norm() < 1.0e-6);
        }
        {
            let tmp0 = &mult_vec(&mat, evec1.as_slice()[..2].try_into().unwrap());
            assert!(tmp0.sub(&evec1.scale(lam1)).norm() < 1.0e-6);
        }
        let evec0 = evec0.scale(1. / lam0.sqrt());
        let evec1 = evec1.scale(1. / lam1.sqrt());
        let sig = [mat[0], mat[1], mat[1], mat[2]];
        let aabb = aabb2(&mat);
        let ndiv = 128;
        let mut sdf_max = f32::MIN;
        for i in 0..ndiv {
            use crate::mat2_col_major::Mat2ColMajor;
            let theta = (i as f32) * 2.0 * std::f32::consts::PI / (ndiv as f32);
            let v = evec0.scale(theta.cos()).add(&evec1.scale(theta.sin()));
            let radrad = sig.mult_vec(&v).dot(&v);
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
    Real: num_traits::Float,
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
    let p_mat: [Real; 6] = [1., 2., 4., 3., 2., 0.];
    let quat0 = crate::quaternion::normalized(&[-3., -2., 0., -1.]);
    let s0_mat = [0.1, 3.0, 1.0];
    let (abc, dabcdt) = wdw_projected_spd_mat3(&p_mat, &quat0, &s0_mat);
    let r0_mat = crate::quaternion::to_mat3_col_major(&quat0);
    let s0_mat = crate::mat3_col_major::from_diagonal(&s0_mat);
    {
        use crate::mat3_col_major::Mat3ColMajor;
        let prs0 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &r0_mat.mult_mat_col_major(&s0_mat),
        );
        let sigma0 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs0,
            &crate::mat2x3_col_major::transpose(&prs0),
        );
        assert!((abc[0] - sigma0[0]).abs() < 1.0e-5);
        assert!((abc[1] - sigma0[1]).abs() < 1.0e-5);
        assert!((abc[2] - sigma0[3]).abs() < 1.0e-5);
    }
    let eps: Real = 1.0e-5;
    for i in 0..3 {
        use crate::mat3_col_major::Mat3ColMajor;
        let mut s1_mat = s0_mat;
        s1_mat[i * 3 + i] += eps;
        let prs1 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &r0_mat.mult_mat_col_major(&s1_mat),
        );
        let sigma1 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs1,
            &crate::mat2x3_col_major::transpose(&prs1),
        );
        let a_si_diff = (sigma1[0] - abc[0]) / eps;
        let b_si_diff = (sigma1[1] - abc[1]) / eps;
        let c_si_diff = (sigma1[3] - abc[2]) / eps;
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
        let r1_mat = crate::mat3_col_major::mult_mat_col_major(&w, &r0_mat);
        let prs1 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &crate::mat3_col_major::mult_mat_col_major(&r1_mat, &s0_mat),
        );
        let sigma1 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs1,
            &crate::mat2x3_col_major::transpose(&prs1),
        );
        let a_ri_diff = (sigma1[0] - abc[0]) / eps;
        let b_ri_diff = (sigma1[1] - abc[1]) / eps;
        let c_ri_diff = (sigma1[3] - abc[2]) / eps;
        let a_ri_ana = dabcdt[0][3 + i];
        let b_ri_ana = dabcdt[1][3 + i];
        let c_ri_ana = dabcdt[2][3 + i];
        assert!((a_ri_diff - a_ri_ana).abs() < 1.0e-3 * (a_ri_ana.abs() + 1.0));
        assert!((c_ri_diff - c_ri_ana).abs() < 1.0e-3 * (c_ri_ana.abs() + 1.0));
        assert!((b_ri_diff - b_ri_ana).abs() < 1.0e-3 * (b_ri_ana.abs() + 1.0));
    }
}

pub fn wdw_inverse<Real, const N: usize>(dabcdt: &[[Real; N]; 3], xyz: &[Real; 3]) -> [[Real; N]; 3]
where
    Real: num_traits::Float,
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
    let p_mat: [Real; 6] = [1., 2., 4., 3., 2., 0.];
    let quat0 = crate::quaternion::normalized(&[-3., -2., 0., -1.]);
    let s0_mat = [0.1, 3.0, 1.0];
    let (abc, dabcdt) = wdw_projected_spd_mat3(&p_mat, &quat0, &s0_mat);
    let xyz = safe_inverse(&abc);
    let dxyzdt = wdw_inverse(&dabcdt, &xyz);
    let r0_mat = crate::quaternion::to_mat3_col_major(&quat0);
    let s0_mat = crate::mat3_col_major::from_diagonal(&s0_mat);
    {
        use crate::mat3_col_major::Mat3ColMajor;
        let prs0 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &r0_mat.mult_mat_col_major(&s0_mat),
        );
        let sigma0 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs0,
            &crate::mat2x3_col_major::transpose(&prs0),
        );
        let sigma0inv = crate::mat2_col_major::try_inverse(&sigma0).unwrap();
        assert!((xyz[0] - sigma0inv[0]).abs() < 1.0e-5);
        assert!((xyz[1] - sigma0inv[1]).abs() < 1.0e-5);
        assert!((xyz[2] - sigma0inv[3]).abs() < 1.0e-5);
    }
    let eps: Real = 1.0e-5;
    for i in 0..3 {
        use crate::mat3_col_major::Mat3ColMajor;
        let mut s1_mat = s0_mat;
        s1_mat[i * 3 + i] += eps;
        let prs1 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &r0_mat.mult_mat_col_major(&s1_mat),
        );
        let sigma1 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs1,
            &crate::mat2x3_col_major::transpose(&prs1),
        );
        let sigma1inv = crate::mat2_col_major::try_inverse(&sigma1).unwrap();
        let x_si_diff = (sigma1inv[0] - xyz[0]) / eps;
        let y_si_diff = (sigma1inv[1] - xyz[1]) / eps;
        let z_si_diff = (sigma1inv[3] - xyz[2]) / eps;
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
        use crate::mat3_col_major::Mat3ColMajor;
        let aa = crate::vec3::basis(i, eps);
        let w = crate::vec3::to_mat3_from_axisangle_vec(&aa);
        let r1_mat = crate::mat3_col_major::mult_mat_col_major(&w, &r0_mat);
        let prs1 = crate::mat2x3_col_major::mult_mat3_col_major(
            &p_mat,
            &r1_mat.mult_mat_col_major(&s0_mat),
        );
        let sigma1 = crate::mat2x3_col_major::mult_mat3x2_col_major(
            &prs1,
            &crate::mat2x3_col_major::transpose(&prs1),
        );
        let sigma1inv = crate::mat2_col_major::try_inverse(&sigma1).unwrap();
        let x_ri_diff = (sigma1inv[0] - xyz[0]) / eps;
        let y_ri_diff = (sigma1inv[1] - xyz[1]) / eps;
        let z_ri_diff = (sigma1inv[3] - xyz[2]) / eps;
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

/// The input can be both column major or row major.
pub fn from_mat2_by_symmetrization<Real>(m: &[Real; 4]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let half = one / (one + one);
    [m[0], (m[1] + m[2]) * half, m[3]]
}
