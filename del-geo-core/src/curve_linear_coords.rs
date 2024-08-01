//! methods for 3x3 matrix

/// covariant to contravariant or contravariant to covariant
pub fn inverse<T>(gd: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    // contravariant basis vectors
    let mut gu = [[T::zero(); 3]; 3];
    crate::vec3::cross_mut(&mut gu[0], &gd[1], &gd[2]);
    let invtmp1 = T::one() / crate::vec3::dot(&gu[0], &gd[0]);
    gu[0][0] = gu[0][0] * invtmp1;
    gu[0][1] = gu[0][1] * invtmp1;
    gu[0][2] = gu[0][2] * invtmp1;
    //
    crate::vec3::cross_mut(&mut gu[1], &gd[2], &gd[0]);
    let invtmp2 = T::one() / crate::vec3::dot(&gu[1], &gd[1]);
    gu[1][0] = gu[1][0] * invtmp2;
    gu[1][1] = gu[1][1] * invtmp2;
    gu[1][2] = gu[1][2] * invtmp2;
    //
    crate::vec3::cross_mut(&mut gu[2], &gd[0], &gd[1]);
    let invtmp3 = T::one() / crate::vec3::dot(&gu[2], &gd[2]);
    gu[2][0] = gu[2][0] * invtmp3;
    gu[2][1] = gu[2][1] * invtmp3;
    gu[2][2] = gu[2][2] * invtmp3;
    gu
}

#[test]
fn hoge() {
    let gds: [[f64; 3]; 3] = [[0., 2., 4.], [3., 5., 4.], [6., 7., 8.]];
    let gus = inverse(&gds);
    for i in 0..3 {
        let gd = &gds[i];
        for j in 0..3 {
            let gu = &gus[j];
            let dot = crate::vec3::dot(gd, gu);
            if i == j {
                assert!((dot - 1.0).abs() < 1.0e-10);
            } else {
                assert!(dot.abs() < 1.0e-10);
            }
        }
    }
}
