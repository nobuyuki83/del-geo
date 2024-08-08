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
    let quat: [f32; 4] = [1., 2., 3., 1.];
    let quat = normalized(&quat);
    let r_mat = to_mat3_col_major(&quat);
    let r_mat_transp = crate::mat3_col_major::transpose(&r_mat);
    {
        let identity0 = crate::mat3_col_major::mult_mat_col_major(&r_mat, &r_mat_transp);
        dbg!(&identity0);
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
    let s_mat = crate::mat3_col_major::from_diagonal(&[1.0, 2.0, 0.3]);
    let rs = crate::mat3_col_major::mult_mat_col_major(&r_mat, &s_mat);
    let a_vec = [0.2, 2.3, 0.1];
    let x_vec = crate::mat3_col_major::mult_vec(&crate::mat3_col_major::transpose(&rs), &a_vec);
    let l = crate::vec3::squared_norm(&x_vec);
    dbg!(&x_vec, l);
}
