pub fn from_identity<T>() -> [T; 9]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, zero, one, zero, zero, zero, one]
}

// above: from method
// --------------------------

pub fn determinant<Real>(u: &[Real; 9]) -> Real
where
    Real: num_traits::Float,
{
    u[0] * u[4] * u[8] + u[3] * u[7] * u[2] + u[6] * u[1] * u[5]
        - u[0] * u[7] * u[5]
        - u[6] * u[4] * u[2]
        - u[3] * u[1] * u[8]
}

pub fn squared_norm<Real>(u: &[Real; 9]) -> Real
where
    Real: num_traits::Float + std::iter::Sum,
{
    u.iter().map(|&v| v * v).sum()
}

pub fn transpose<Real>(m: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

pub fn mult_mat_row_major<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float + std::ops::AddAssign,
{
    let mut r = [Real::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j]
            }
        }
    }
    r
}

pub fn sub<Real>(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9]
where
    Real: num_traits::Float,
{
    [
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2],
        a[3] - b[3],
        a[4] - b[4],
        a[5] - b[5],
        a[6] - b[6],
        a[7] - b[7],
        a[8] - b[8],
    ]
}

fn sort_eigen(g: &mut [f64; 3], v: &mut [f64; 9]) {
    if g[1] > g[0] {
        // swap 01
        let mut t = g[0];
        g[0] = g[1];
        g[1] = t;
        t = v[0];
        v[0] = v[1];
        v[1] = t;
        t = v[3];
        v[3] = v[4];
        v[4] = t;
        t = v[6];
        v[6] = v[7];
        v[7] = t;
    }
    if g[2] > g[1] {
        let mut t = g[1];
        g[1] = g[2];
        g[2] = t;
        t = v[1];
        v[1] = v[2];
        v[2] = t;
        t = v[4];
        v[4] = v[5];
        v[5] = t;
        t = v[7];
        v[7] = v[8];
        v[8] = t;
    }
    if g[1] > g[0] {
        // swap 01
        let mut t = g[0];
        g[0] = g[1];
        g[1] = t;
        t = v[0];
        v[0] = v[1];
        v[1] = t;
        t = v[3];
        v[3] = v[4];
        v[4] = t;
        t = v[6];
        v[6] = v[7];
        v[7] = t;
    }
}

// m = UGV^T
pub fn svd(m: &[f64; 9], nitr: usize) {
    // M^TM = VGGV^T
    let mtm = [
        m[0] * m[0] + m[3] * m[3] + m[6] * m[6],
        m[1] * m[1] + m[4] * m[4] + m[7] * m[7],
        m[2] * m[2] + m[5] * m[5] + m[8] * m[8],
        m[1] * m[2] + m[4] * m[5] + m[7] * m[8],
        m[2] * m[0] + m[5] * m[3] + m[8] * m[6],
        m[0] * m[1] + m[3] * m[4] + m[6] * m[7],
    ];
    let Some((mut v, mut lv)) = crate::mat3_sym::eigen_decomp(mtm, nitr) else {
        todo!()
    };
    sort_eigen(&mut lv, &mut v);
    if lv[0] < 0. {
        lv[0] = 0.0;
    }
    if lv[1] < 0. {
        lv[1] = 0.0;
    }
    if lv[2] < 0. {
        lv[2] = 0.0;
    }
    let mut g = [lv[0].sqrt(), lv[1].sqrt(), lv[2].sqrt()];

    let mut u0 = [
        m[0] * v[0] + m[1] * v[3] + m[2] * v[6],
        m[3] * v[0] + m[4] * v[3] + m[5] * v[6],
        m[6] * v[0] + m[7] * v[3] + m[8] * v[6],
    ];
    let mut u1 = [
        m[0] * v[1] + m[1] * v[4] + m[2] * v[7],
        m[3] * v[1] + m[4] * v[4] + m[5] * v[7],
        m[6] * v[1] + m[7] * v[4] + m[8] * v[7],
    ];
    let mut u2 = [
        m[0] * v[2] + m[1] * v[5] + m[2] * v[8],
        m[3] * v[2] + m[4] * v[5] + m[5] * v[8],
        m[6] * v[2] + m[7] * v[5] + m[8] * v[8],
    ];

    let mut u = [0f64; 9];
    if determinant(&v) < 0. {
        // making right hand coordinate
        v[0 * 3 + 2] *= -1.;
        v[1 * 3 + 2] *= -1.;
        v[2 * 3 + 2] *= -1.;
        g[2] *= -1.;
        u[0 * 3 + 2] *= -1.;
        u[1 * 3 + 2] *= -1.;
        u[2 * 3 + 2] *= -1.;
    }

    let sql0 = crate::vec3::squared_norm(&u0);
    if sql0 > 1.0e-20 {
        crate::vec3::normalize(&mut u0);
        let sql1 = crate::vec3::squared_norm(&u1);
        if sql1 < 1.0e-20 {
            u1[0] = 1.0 - u0[0].abs();
            u1[1] = 1.0 - u0[1].abs();
            u1[2] = 1.0 - u0[2].abs();
        } else {
            crate::vec3::normalize(&mut u1);
        }
        let d01 = u0[0] * u1[0] + u0[1] * u1[1] + u0[2] * u1[2];
        u1[0] -= d01 * u0[0];
        u1[1] -= d01 * u0[1];
        u1[2] -= d01 * u0[2];
        crate::vec3::normalize(&mut u1);
        let s2 = crate::vec3::cross(&u0, &u1);
        let d22 = u2[0] * s2[0] + u2[1] * s2[1] + u2[2] * s2[2];
        u2[0] = s2[0];
        u2[1] = s2[1];
        u2[2] = s2[2];
        if d22 < 0. {
            g[2] *= -1.;
        }
    } else {
        u0[0] = 1.;
        u0[1] = 0.;
        u0[2] = 0.;
        u1[0] = 0.;
        u1[1] = 1.;
        u1[2] = 0.;
        u2[0] = 0.;
        u2[1] = 0.;
        u2[2] = 1.;
    }
    u[0] = u0[0];
    u[1] = u1[0];
    u[2] = u2[0];
    u[3] = u0[1];
    u[4] = u1[1];
    u[5] = u2[1];
    u[6] = u0[2];
    u[7] = u1[2];
    u[8] = u2[2];
}
