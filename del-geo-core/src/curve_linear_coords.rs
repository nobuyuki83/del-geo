//! methods for 3x3 matrix

use crate::vec3::Vec3;

/// covariant to contravariant or contravariant to covariant
pub fn inverse<T>(gd: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float + std::ops::MulAssign,
{
    // contravariant basis vectors
    let mut gu = [[T::zero(); 3]; 3];
    gu[0].cross_mut(&gd[1], &gd[2]);
    let invtmp1 = T::one() / gu[0].dot(&gd[0]);
    gu[0] = gu[0].map(|x| x * invtmp1);
    //
    gu[1].cross_mut(&gd[2], &gd[0]);
    let invtmp2 = T::one() / gu[1].dot(&gd[1]);
    gu[1] = gu[1].map(|x| x * invtmp2);
    //
    gu[2].cross_mut(&gd[0], &gd[1]);
    let invtmp3 = T::one() / gu[2].dot(&gd[2]);
    gu[2] = gu[2].map(|x| x * invtmp3);
    gu
}

#[test]
fn hoge() {
    let gds = [[0f64, 2., 4.], [3., 5., 4.], [6., 7., 8.]];
    let gus = inverse(&gds);
    for (i, gd) in gds.iter().enumerate() {
        for (j, gu) in gus.iter().enumerate() {
            let dot = gd.dot(gu);
            if i == j {
                assert!((dot - 1.0).abs() < 1.0e-10);
            } else {
                assert!(dot.abs() < 1.0e-10);
            }
        }
    }
}
