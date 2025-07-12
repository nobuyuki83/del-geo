//! methods for 2D triangle (parameterized by three corner points)

pub fn area<T>(v1: &nalgebra::Vector2<T>, v2: &nalgebra::Vector2<T>, v3: &nalgebra::Vector2<T>) -> T
where
    T: nalgebra::RealField + Copy,
{
    ((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1])) / (T::one() + T::one())
}

pub fn circumcenter<T>(
    p0: &nalgebra::Vector2<T>,
    p1: &nalgebra::Vector2<T>,
    p2: &nalgebra::Vector2<T>,
) -> nalgebra::Vector2<T>
where
    T: nalgebra::RealField + Copy + std::fmt::Debug,
{
    let a0 = (p1 - p2).norm_squared();
    let a1 = (p2 - p0).norm_squared();
    let a2 = (p0 - p1).norm_squared();

    let b0: T = a0 * (a1 + a2 - a0);
    let b1: T = a1 * (a0 + a2 - a1);
    let b2: T = a2 * (a0 + a1 - a2);
    let sum = T::one() / (b0 + b1 + b2);

    let c0 = b0 * sum;
    let c1 = b1 * sum;
    let c2 = b2 * sum;

    p0.scale(c0) + p1.scale(c1) + p2.scale(c2)
}

#[test]
fn test_circumcenter() {
    type VEC = nalgebra::Vector2<f64>;
    let p0 = [VEC::new(1.3, 2.1), VEC::new(3.2, 2.1), VEC::new(1.5, 2.5)];
    let cc0 = circumcenter(&p0[0], &p0[1], &p0[2]);
    let d0 = (cc0 - p0[0]).norm_squared();
    let d1 = (cc0 - p0[1]).norm_squared();
    let d2 = (cc0 - p0[2]).norm_squared();
    assert!((d0 - d1).abs() < d0 * 1.0e-10);
    assert!((d0 - d2).abs() < d0 * 1.0e-10);
}

pub fn wdw_circumcenter<T>(
    p0: &nalgebra::Vector2<T>,
    p1: &nalgebra::Vector2<T>,
    p2: &nalgebra::Vector2<T>,
) -> (nalgebra::Vector2<T>, [nalgebra::Matrix2<T>; 3])
where
    T: nalgebra::RealField + Copy + std::fmt::Debug + num_traits::Float,
{
    let a0 = (p1 - p2).norm_squared();
    let a1 = (p2 - p0).norm_squared();
    let a2 = (p0 - p1).norm_squared();
    //
    let b0: T = a0 * (a1 + a2 - a0);
    let b1: T = a1 * (a2 + a0 - a1);
    let b2: T = a2 * (a0 + a1 - a2);
    //
    let sum = b0 + b1 + b2;
    let sum_inv = T::one() / sum;
    //
    let c0 = b0 * sum_inv;
    let c1 = b1 * sum_inv;
    let c2 = b2 * sum_inv;
    let cc = p0.scale(c0) + p1.scale(c1) + p2.scale(c2);
    // -----------------
    let two = T::one() + T::one();
    let db0 = [
        (p0 - p2 + p0 - p1).scale(two * a0),
        (p1 - p0 + p2 - p1).scale(two * a0) + (p1 - p2).scale(two * (a1 + a2 - a0)),
        (p2 - p0 + p1 - p2).scale(two * a0) + (p2 - p1).scale(two * (a1 + a2 - a0)),
    ];
    let db1 = [
        (p0 - p1 + p2 - p0).scale(two * a1) + (p0 - p2).scale(two * (a2 + a0 - a1)),
        (p1 - p0 + p1 - p2).scale(two * a1),
        (p2 - p1 + p0 - p2).scale(two * a1) + (p2 - p0).scale(two * (a2 + a0 - a1)),
    ];
    let db2 = [
        (p0 - p2 + p1 - p0).scale(two * a2) + (p0 - p1).scale(two * (a0 + a1 - a2)),
        (p1 - p2 + p0 - p1).scale(two * a2) + (p1 - p0).scale(two * (a0 + a1 - a2)),
        (p2 - p1 + p2 - p0).scale(two * a2),
    ];
    let tmp = -T::one() / (sum * sum);
    let dsum_inv = [
        (db0[0] + db1[0] + db2[0]) * tmp,
        (db0[1] + db1[1] + db2[1]) * tmp,
        (db0[2] + db1[2] + db2[2]) * tmp,
    ];
    //
    let dcc = [
        nalgebra::Matrix2::<T>::identity() * c0
            + p0 * (db0[0].scale(sum_inv) + dsum_inv[0].scale(b0)).transpose()
            + p1 * (db1[0].scale(sum_inv) + dsum_inv[0].scale(b1)).transpose()
            + p2 * (db2[0].scale(sum_inv) + dsum_inv[0].scale(b2)).transpose(),
        nalgebra::Matrix2::<T>::identity() * c1
            + p0 * (db0[1].scale(sum_inv) + dsum_inv[1].scale(b0)).transpose()
            + p1 * (db1[1].scale(sum_inv) + dsum_inv[1].scale(b1)).transpose()
            + p2 * (db2[1].scale(sum_inv) + dsum_inv[1].scale(b2)).transpose(),
        nalgebra::Matrix2::<T>::identity() * c2
            + p0 * (db0[2].scale(sum_inv) + dsum_inv[2].scale(b0)).transpose()
            + p1 * (db1[2].scale(sum_inv) + dsum_inv[2].scale(b1)).transpose()
            + p2 * (db2[2].scale(sum_inv) + dsum_inv[2].scale(b2)).transpose(),
    ];

    (cc, dcc)
}

#[test]
fn test_dw_circumcenter() {
    type VEC = nalgebra::Vector2<f64>;
    let p0 = [VEC::new(0.1, 0.2), VEC::new(1.3, 0.2), VEC::new(0.3, 1.5)];
    let (cc0, dcc0) = wdw_circumcenter(&p0[0], &p0[1], &p0[2]);
    let eps = 1.0e-4;
    for i_node in 0..3 {
        for i_dim in 0..2 {
            let p1 = {
                let mut p1 = p0;
                p1[i_node][i_dim] += eps;
                p1
            };
            let (cc1, _dcc1) = wdw_circumcenter(&p1[0], &p1[1], &p1[2]);
            let dcc_num = (cc1 - cc0) / eps;
            let mut b = nalgebra::Vector2::<f64>::zeros();
            b[i_dim] = 1.0;
            let dcc_ana = dcc0[i_node] * b;
            let diff = (dcc_num - dcc_ana).norm();
            assert!(diff < 1.0e-4);
            println!("{i_node}, {i_dim} --> {dcc_num:?}, {dcc_ana:?}, {diff:?}");
        }
    }
}
